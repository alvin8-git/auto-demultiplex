#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import gzip
import argparse
from collections import defaultdict, Counter

prog_version = '1.1.0'
prog_date = '2024-01-01'

DEFAULT_BARCODE_FILES = [
    'BarcodeV2.1.txt',
    'BarcodeV3.0.txt',
]

def rev_comp(s):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(comp.get(c, c) for c in s[::-1])

def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        return max(len(s1), len(s2))
    return sum(c1 != c2 and c1 != 'N' and c2 != 'N' for c1, c2 in zip(s1, s2))

def get_segment_lengths(barcode_seq, seg_len=10):
    """Determine individual barcode segment lengths from a (possibly N-padded) sequence.

    Rules (all barcodes in a run use fixed segment sizes of seg_len bp):
    - All non-N: split at seg_len boundaries (e.g. 10 → [10], 20 → [10,10], 16 → [10,6])
    - Non-N prefix then all-N suffix: first segment = len(non-N prefix), second = rest
    - All-N prefix then non-N suffix: first segment = len(N prefix), second = rest
    """
    first_n = next((i for i, c in enumerate(barcode_seq) if c == 'N'), None)
    if first_n is None:
        # All non-N: split at seg_len boundaries
        segs, pos = [], 0
        while pos < len(barcode_seq):
            segs.append(min(seg_len, len(barcode_seq) - pos))
            pos += seg_len
        return segs
    if first_n == 0:
        # Starts with N: first segment = N-prefix, second = rest
        first_non_n = next((i for i, c in enumerate(barcode_seq) if c != 'N'), len(barcode_seq))
        return [first_non_n, len(barcode_seq) - first_non_n]
    # Starts with non-N, then N: first segment = non-N prefix, second = rest
    return [first_n, len(barcode_seq) - first_n]

def apply_segment_rc(read_bc, seg_lengths):
    """RC each barcode segment independently and concatenate."""
    result, pos = [], 0
    for seg_len in seg_lengths:
        result.append(rev_comp(read_bc[pos:pos + seg_len]))
        pos += seg_len
    return ''.join(result)

def validate_input_file(path, name):
    if not os.path.exists(path):
        print(f'ERROR: {name} file not found: {path}', file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(path):
        print(f'ERROR: {name} is not a file: {path}', file=sys.stderr)
        sys.exit(1)
    if not os.access(path, os.R_OK):
        print(f'ERROR: {name} is not readable: {path}', file=sys.stderr)
        sys.exit(1)
    return True

def load_barcode_file(barcode_path):
    barcodes = {}
    if not os.path.exists(barcode_path):
        return barcodes
    with open(barcode_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                name = parts[0]
                seq = parts[1].upper()
                barcodes[name] = seq
    return barcodes

def load_all_known_barcodes():
    all_barcodes = {}
    for bc_file in DEFAULT_BARCODE_FILES:
        if os.path.exists(bc_file):
            barcodes = load_barcode_file(bc_file)
            all_barcodes.update(barcodes)
    return all_barcodes

def read_fastq_records(fq_path, sample_size=None, progress=False):
    records = []
    total_reads = 0
    
    f = gzip.open(fq_path, 'rb') if fq_path.endswith('.gz') else open(fq_path, 'rb')
    
    while True:
        lines = []
        for _ in range(4):
            l = f.readline()
            if not l:
                break
            lines.append(l)
        
        if not lines:
            break
        
        seq = lines[1]
        if isinstance(seq, bytes):
            seq = seq.decode()
        records.append(seq.strip())
        
        total_reads += 1
        
        if sample_size and total_reads >= sample_size:
            break
        
        if progress and total_reads % 10000 == 0:
            print(f'  Read {total_reads} records...', end='\r')
    
    f.close()
    
    if progress:
        print(f'  Read {total_reads} records complete    ')
    
    return records

def detect_barcode_config(r1_path, r2_path=None, sample_size=50000):
    print('Detecting barcode configuration...')
    
    known_barcodes = load_all_known_barcodes()
    print(f'Loaded {len(known_barcodes)} known barcodes')
    
    if r2_path:
        return detect_pe_config(r1_path, r2_path, known_barcodes, sample_size)
    else:
        return detect_se_config(r1_path, known_barcodes, sample_size)

def find_matching_barcode(seq, known_barcodes, max_mismatch=1):
    for bc, count in Counter([seq, rev_comp(seq)]).items():
        for name, known_seq in known_barcodes.items():
            dist = hamming_distance(bc, known_seq)
            if dist <= max_mismatch:
                return name, known_seq, bc == rev_comp(known_seq), dist
    return None, None, None, None

def detect_pe_config(r1_path, r2_path, known_barcodes, sample_size):
    print(f'Scanning {r2_path} for barcode patterns...')
    
    r2_seqs = read_fastq_records(r2_path, sample_size, progress=True)
    
    if not r2_seqs:
        print('ERROR: No reads found in R2 file', file=sys.stderr)
        sys.exit(1)
    
    read_len = len(r2_seqs[0])
    barcodes_counter = Counter()
    
    for seq in r2_seqs:
        for bc_len in [10, 20]:
            if len(seq) >= bc_len:
                barcodes_counter[seq[-bc_len:]] += 1
    
    matched = {}
    for bc, count in barcodes_counter.most_common(200):
        name, known_seq, is_rc, dist = find_matching_barcode(bc, known_barcodes)
        if name:
            matched[name] = {
                'sequence': known_seq,
                'use_rc': is_rc,
                'mismatch': dist,
                'count': count
            }
    
    if matched:
        bc_len = len(next(iter(matched.values()))['sequence'])
        return {
            'mode': 'PE',
            'barcodes': matched,
            'barcode_length': bc_len,
            'read_len_before_barcode': read_len - bc_len,
            'use_reverse': True,
        }
    
    print('No known barcodes found, extracting from data...')
    return extract_unknown_pe(r2_seqs, read_len)

def extract_unknown_pe(r2_seqs, read_len):
    barcodes_10 = Counter()
    
    for seq in r2_seqs:
        if len(seq) >= 10:
            barcodes_10[seq[-10:]] += 1
    
    filtered = {bc: count for bc, count in barcodes_10.items() 
                if count >= 10 and 'N' not in bc}
    
    extracted = {}
    for i, (bc, count) in enumerate(sorted(filtered.items(), key=lambda x: -x[1])[:100]):
        name = f'barcode{i+1}'
        extracted[name] = {
            'sequence': bc,
            'use_rc': False,
            'mismatch': 0,
            'count': count
        }
    
    bc_len = 10
    return {
        'mode': 'PE',
        'barcodes': extracted,
        'barcode_length': bc_len,
        'read_len_before_barcode': read_len - bc_len,
        'use_reverse': False,
    }

def detect_se_config(fq_path, known_barcodes, sample_size):
    print(f'Scanning {fq_path} for barcode patterns...')
    
    fq_seqs = read_fastq_records(fq_path, sample_size, progress=True)
    
    if not fq_seqs:
        print('ERROR: No reads found in FASTQ file', file=sys.stderr)
        sys.exit(1)
    
    read_len = len(fq_seqs[0])
    barcodes_counter = Counter()
    
    for seq in fq_seqs:
        for bc_len in [10, 20]:
            if len(seq) >= bc_len:
                barcodes_counter[seq[-bc_len:]] += 1

    matched = {}
    for bc, count in barcodes_counter.most_common(200):
        name, known_seq, is_rc, dist = find_matching_barcode(bc, known_barcodes)
        if name:
            matched[name] = {
                'sequence': known_seq,
                'use_rc': is_rc,
                'mismatch': dist,
                'count': count
            }

    if matched:
        bc_len = len(next(iter(matched.values()))['sequence'])
        return {
            'mode': 'SE',
            'barcodes': matched,
            'barcode_length': bc_len,
            'read_len_before_barcode': read_len - bc_len,
            'use_reverse': False,
        }
    
    return extract_unknown_se(fq_seqs, read_len)

def extract_unknown_se(fq_seqs, read_len):
    barcodes_10 = Counter()
    
    for seq in fq_seqs:
        if len(seq) >= 10:
            barcodes_10[seq[-10:]] += 1
    
    filtered = {bc: count for bc, count in barcodes_10.items() 
                if count >= 10 and 'N' not in bc}
    
    extracted = {}
    for i, (bc, count) in enumerate(sorted(filtered.items(), key=lambda x: -x[1])[:100]):
        name = f'barcode{i+1}'
        extracted[name] = {
            'sequence': bc,
            'use_rc': False,
            'mismatch': 0,
            'count': count
        }
    
    bc_len = 10
    return {
        'mode': 'SE',
        'barcodes': extracted,
        'barcode_length': bc_len,
        'read_len_before_barcode': read_len - bc_len,
        'use_reverse': False,
    }

def demultiplex_pe(r1_path, r2_path, barcodes, output_dir, config, max_mismatch=2):
    print(f'Demultiplexing {r1_path} + {r2_path}...')

    os.makedirs(output_dir, exist_ok=True)

    bc_len = config['barcode_length']
    use_rc = config['use_reverse']

    # Determine per-segment RC boundaries from first barcode sequence
    first_seq = next(iter(barcodes.values()))['sequence']
    seg_lengths = get_segment_lengths(first_seq)

    # Exact-match lookup (only for barcodes without N wildcards)
    barcode_lookup = {}
    for name, info in barcodes.items():
        seq = info['sequence']
        if 'N' not in seq:
            barcode_lookup[seq] = name

    barcode_files_r1 = {}
    barcode_files_r2 = {}

    for name in barcodes:
        fq1_out = gzip.open(os.path.join(output_dir, f'{name}_1.fq.gz'), 'wb')
        fq2_out = gzip.open(os.path.join(output_dir, f'{name}_2.fq.gz'), 'wb')
        barcode_files_r1[name] = fq1_out
        barcode_files_r2[name] = fq2_out

    ambiguous_r1 = gzip.open(os.path.join(output_dir, 'ambiguous_1.fq.gz'), 'wb')
    ambiguous_r2 = gzip.open(os.path.join(output_dir, 'ambiguous_2.fq.gz'), 'wb')

    f1 = gzip.open(r1_path, 'rb') if r1_path.endswith('.gz') else open(r1_path, 'rb')
    f2 = gzip.open(r2_path, 'rb') if r2_path.endswith('.gz') else open(r2_path, 'rb')

    stats = defaultdict(int)
    mismatch_stats = defaultdict(int)
    seq_stats = Counter()

    total_reads = 0

    while True:
        lines1 = []
        lines2 = []

        for _ in range(4):
            l1 = f1.readline()
            l2 = f2.readline()
            if not l1 or not l2:
                break
            lines1.append(l1)
            lines2.append(l2)

        if not lines1 or not lines1[0]:
            break

        seq2 = lines2[1]
        if isinstance(seq2, bytes):
            seq2 = seq2.decode()
        seq2 = seq2.strip()

        raw_bc = seq2[-bc_len:]
        bc = apply_segment_rc(raw_bc, seg_lengths) if use_rc else raw_bc

        matched_name = barcode_lookup.get(bc)

        if not matched_name:
            best_name, best_dist = None, max_mismatch + 1
            for name, info in barcodes.items():
                dist = hamming_distance(bc, info['sequence'])
                if dist < best_dist:
                    best_dist = dist
                    best_name = name
            if best_name is not None and best_dist <= max_mismatch:
                matched_name = best_name
                mismatch_stats[best_name] += 1

        if matched_name:
            barcode_files_r1[matched_name].writelines(lines1)
            barcode_files_r2[matched_name].writelines(lines2)
            stats[matched_name] += 1
            seq_stats[matched_name] += 1
        else:
            ambiguous_r1.writelines(lines1)
            ambiguous_r2.writelines(lines2)
            stats['ambiguous'] += 1

        stats['total'] += 1
        total_reads += 1

        if total_reads % 50000 == 0:
            print(f'  Processed {total_reads} reads...', end='\r')

    f1.close()
    f2.close()

    for f in barcode_files_r1.values():
        f.close()
    for f in barcode_files_r2.values():
        f.close()
    ambiguous_r1.close()
    ambiguous_r2.close()

    print(f'  Processed {total_reads} reads complete    ')

    return stats, seq_stats, mismatch_stats

def demultiplex_se(fq_path, barcodes, output_dir, config, max_mismatch=2):
    print(f'Demultiplexing {fq_path}...')

    os.makedirs(output_dir, exist_ok=True)

    bc_len = config['barcode_length']
    use_rc = config['use_reverse']

    first_seq = next(iter(barcodes.values()))['sequence']
    seg_lengths = get_segment_lengths(first_seq)

    barcode_lookup = {}
    for name, info in barcodes.items():
        seq = info['sequence']
        if 'N' not in seq:
            barcode_lookup[seq] = name

    barcode_files = {}
    for name in barcodes:
        fq_out = gzip.open(os.path.join(output_dir, f'{name}.fq.gz'), 'wb')
        barcode_files[name] = fq_out

    ambiguous_file = gzip.open(os.path.join(output_dir, 'ambiguous.fq.gz'), 'wb')

    f = gzip.open(fq_path, 'rb') if fq_path.endswith('.gz') else open(fq_path, 'rb')

    stats = defaultdict(int)
    mismatch_stats = defaultdict(int)
    seq_stats = Counter()

    total_reads = 0

    while True:
        lines = []
        for _ in range(4):
            l = f.readline()
            if not l:
                break
            lines.append(l)

        if not lines or not lines[0]:
            break

        read_seq = lines[1]
        if isinstance(read_seq, bytes):
            read_seq = read_seq.decode()
        read_seq = read_seq.strip()

        raw_bc = read_seq[-bc_len:]
        bc = apply_segment_rc(raw_bc, seg_lengths) if use_rc else raw_bc

        matched_name = barcode_lookup.get(bc)

        if not matched_name:
            best_name, best_dist = None, max_mismatch + 1
            for name, info in barcodes.items():
                dist = hamming_distance(bc, info['sequence'])
                if dist < best_dist:
                    best_dist = dist
                    best_name = name
            if best_name is not None and best_dist <= max_mismatch:
                matched_name = best_name
                mismatch_stats[best_name] += 1

        if matched_name:
            barcode_files[matched_name].writelines(lines)
            stats[matched_name] += 1
            seq_stats[matched_name] += 1
        else:
            ambiguous_file.writelines(lines)
            stats['ambiguous'] += 1

        stats['total'] += 1
        total_reads += 1

        if total_reads % 50000 == 0:
            print(f'  Processed {total_reads} reads...', end='\r')

    f.close()

    for f in barcode_files.values():
        f.close()
    ambiguous_file.close()

    print(f'  Processed {total_reads} reads complete    ')

    return stats, seq_stats, mismatch_stats

def write_statistics(barcodes, stats, seq_stats, mismatch_stats, output_dir):
    total = stats['total']
    demux_total = sum(v for k, v in stats.items() if k != 'total' and k != 'ambiguous')
    
    with open(os.path.join(output_dir, 'BarcodeStat.txt'), 'w') as f:
        f.write('BarcodeID\tPerfectMatch\tMismatchMatch\tTotalReads\tPercent\n')
        for name in sorted(barcodes.keys()):
            count = stats.get(name, 0)
            mismatch = mismatch_stats.get(name, 0)
            pct = 100.0 * count / total if total > 0 else 0
            f.write(f'{name}\t{count - mismatch}\t{mismatch}\t{count}\t{pct:.2f}\n')
        
        f.write(f'Total\t{demux_total - sum(mismatch_stats.values())}\t{sum(mismatch_stats.values())}\t{demux_total}\t{100.0*demux_total/total:.2f}\n')
    
    with open(os.path.join(output_dir, 'SequenceStat.txt'), 'w') as f:
        f.write('BarcodeSequence\tBarcodeID\tCount\tPercent\n')
        for name, info in sorted(barcodes.items()):
            seq = info['sequence']
            count = seq_stats.get(name, 0)
            pct = 100.0 * count / total if total > 0 else 0
            f.write(f'{seq}\t{name}\t{count}\t{pct:.2f}\n')

def generate_barcode_file(barcodes, output_path):
    with open(output_path, 'w') as f:
        for name, info in sorted(barcodes.items()):
            f.write(f'{name}\t{info["sequence"]}\n')

def main():
    usage = '''%(prog)s v{version} - Auto-detecting demultiplexing tool for MGI sequencing data

This tool mimics splitBarcode behavior with auto-detection:
  - Barcodes at end of reads (SE) or end of R2 (PE)
  - Reverse complement matching for PE mode
  - Auto-detect barcode length and position
  - Supports mismatch tolerance

Usage:
  %(prog)s -1 <fastq> -o <output_dir>                   # Auto-detect SE
  %(prog)s -1 <r1.fq> -2 <r2.fq> -o <output_dir>       # Auto-detect PE
  %(prog)s -1 <fastq> -B <barcode_file> -o <output_dir> # Use provided barcodes
'''
    
    parser = argparse.ArgumentParser(usage=usage.format(version=prog_version))
    parser.add_argument('-v', '--version', action='version', version=f'{prog_version} ({prog_date})')
    parser.add_argument('-1', '--fastq1', required=True, help='Input FASTQ file (R1 for PE)')
    parser.add_argument('-2', '--fastq2', help='Input FASTQ file R2 (for paired-end)')
    parser.add_argument('-B', '--barcode', help='Barcode file (tsv format: name\tsequence)')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-r', '--reverse', action='store_true', help='Apply reverse complement')
    parser.add_argument('-m', '--mismatch', type=int, default=2, help='Allowed mismatches (default: 2, per MGI recommendation for 10bp barcodes)')
    parser.add_argument('--sample_size', type=int, default=50000, help='Number of reads to sample')
    
    args = parser.parse_args()
    
    validate_input_file(args.fastq1, 'Input FASTQ')
    if args.fastq2:
        validate_input_file(args.fastq2, 'R2 FASTQ')
    if args.barcode:
        validate_input_file(args.barcode, 'Barcode')
    
    output_dir = args.output
    
    if args.barcode:
        raw_barcodes = load_barcode_file(args.barcode)
        if not raw_barcodes:
            print(f'ERROR: No barcodes found in {args.barcode}', file=sys.stderr)
            sys.exit(1)
        barcodes = {name: {'sequence': seq, 'use_rc': args.reverse, 'mismatch': 0, 'count': 0} 
                    for name, seq in raw_barcodes.items()}
        print(f'Loaded {len(barcodes)} barcodes from {args.barcode}')
        
        r2_path = args.fastq2 if args.fastq2 else args.fastq1
        sample = read_fastq_records(r2_path, 100)
        read_len = len(sample[0]) if sample else 100
        
        bc_len = len(next(iter(barcodes.values()))['sequence'])
        config = {
            'mode': 'PE' if args.fastq2 else 'SE',
            'barcode_length': bc_len,
            'use_reverse': args.reverse,
            'read_len_before_barcode': read_len - bc_len,
        }
    else:
        config = detect_barcode_config(args.fastq1, args.fastq2, args.sample_size)
        barcodes = config['barcodes']
        
        print(f'\nDetected configuration:')
        print(f'  Mode: {config["mode"]}')
        print(f'  Barcode length: {config["barcode_length"]}')
        print(f'  Read length before barcode: {config["read_len_before_barcode"]}')
        print(f'  Use reverse complement: {config["use_reverse"]}')
        print(f'  Barcodes found: {len(barcodes)}')
    
    if not barcodes:
        print('ERROR: No barcodes found!', file=sys.stderr)
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)
    barcode_file = os.path.join(output_dir, 'used_barcodes.txt')
    generate_barcode_file(barcodes, barcode_file)
    print(f'Generated barcode file: {barcode_file}')
    
    if config['mode'] == 'PE':
        stats, seq_stats, mismatch_stats = demultiplex_pe(
            args.fastq1, args.fastq2, barcodes, output_dir, config, args.mismatch
        )
    else:
        stats, seq_stats, mismatch_stats = demultiplex_se(
            args.fastq1, barcodes, output_dir, config, args.mismatch
        )
    
    write_statistics(barcodes, stats, seq_stats, mismatch_stats, output_dir)
    
    print(f'\nDemultiplexing complete. Total reads: {stats["total"]}')
    demux_total = sum(v for k, v in stats.items() if k != 'total' and k != 'ambiguous')
    print(f'Demultiplexed: {demux_total} ({100.0*demux_total/stats["total"]:.1f}%)')
    print(f'Ambiguous: {stats.get("ambiguous", 0)}')
    
    for name in sorted(barcodes.keys()):
        count = stats.get(name, 0)
        pct = 100.0 * count / stats['total'] if stats['total'] > 0 else 0
        mismatch = mismatch_stats.get(name, 0)
        print(f'  {name}: {count} ({pct:.1f}%)', end='')
        if mismatch > 0:
            print(f' [with {mismatch} mismatches]')
        else:
            print()
    
    print(f'\nOutput written to: {output_dir}')
    print('Done!')

if __name__ == '__main__':
    main()
