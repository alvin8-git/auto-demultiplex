#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import subprocess
import shutil
import gzip

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

def run_cmd(cmd, check=True, cwd=None):
    print(f'Running: {" ".join(cmd)}')
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd or SCRIPT_DIR)
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)
    if check and result.returncode != 0:
        print(f'FAILED with code {result.returncode}')
        return False
    return True

def count_reads(fq_path):
    if not os.path.exists(fq_path):
        return 0
    f = gzip.open(fq_path, 'rb') if fq_path.endswith('.gz') else open(fq_path, 'rb')
    count = sum(1 for _ in f) // 4
    f.close()
    return count

def total_demuxed(output_dir, barcode_names):
    """Count total reads across all named barcode output files (R1)."""
    total = 0
    for name in barcode_names:
        total += count_reads(os.path.join(output_dir, f'{name}_1.fq.gz'))
    return total

def read_barcode_stat(output_dir):
    stat_path = os.path.join(output_dir, 'BarcodeStat.txt')
    if not os.path.exists(stat_path):
        return {}
    result = {}
    with open(stat_path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4 and parts[0] != 'Total':
                result[parts[0]] = {
                    'perfect': int(parts[1]),
                    'mismatch': int(parts[2]),
                    'total': int(parts[3]),
                }
    return result

def read_sequence_stat(output_dir):
    stat_path = os.path.join(output_dir, 'SequenceStat.txt')
    if not os.path.exists(stat_path):
        return {}
    result = {}
    with open(stat_path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                result[parts[1]] = int(parts[2])  # barcode_id -> count
    return result


# ─── UNIT TESTS ──────────────────────────────────────────────────────────────

def test_unit_hamming_n_wildcard():
    """N in either sequence should count as 0 mismatches."""
    from auto_demultiplex import hamming_distance
    assert hamming_distance('AAAAAAAAAA', 'NNNNNNNNNN') == 0, \
        'N in reference should match any base'
    assert hamming_distance('NNNNNNNNNN', 'AAAAAAAAAA') == 0, \
        'N in read should match any base'
    assert hamming_distance('TGGGAACGCTAAAAAAAAAA', 'TGGGAACGCTNNNNNNNNNN') == 0, \
        'N in second half of reference should not count as mismatch'
    assert hamming_distance('AAAAAAAAAACCGAGACGTT', 'NNNNNNNNNNCCGAGACGTT') == 0, \
        'N in first half of reference should not count as mismatch'
    assert hamming_distance('TGGGAACGCTAAAAAAAAAT', 'TGGGAACGCTNNNNNNNNNA') == 1, \
        'Single real mismatch (last pos T vs A) should still count even with surrounding Ns'
    print('PASSED')
    return True

def test_unit_hamming_normal():
    """Normal hamming distance (no N) should still work correctly."""
    from auto_demultiplex import hamming_distance
    assert hamming_distance('AAAA', 'AAAA') == 0
    assert hamming_distance('AAAA', 'AAAT') == 1
    assert hamming_distance('AAAA', 'TTTT') == 4
    # Unequal lengths still return max
    assert hamming_distance('AAAA', 'AAAAA') == 5
    print('PASSED')
    return True

def test_unit_get_segment_lengths_single():
    """Single 10bp barcode (all non-N) → one segment of 10."""
    from auto_demultiplex import get_segment_lengths
    assert get_segment_lengths('TGGGAACGCT') == [10], \
        'Single 10bp non-N barcode should give [10]'
    print('PASSED')
    return True

def test_unit_get_segment_lengths_dual_10_10():
    """Dual 10+10bp (all non-N 20bp) → two segments of 10."""
    from auto_demultiplex import get_segment_lengths
    assert get_segment_lengths('TGGGAACGCTCCGAGACGTT') == [10, 10], \
        'Dual 20bp non-N barcode should give [10, 10]'
    print('PASSED')
    return True

def test_unit_get_segment_lengths_dual_10_6():
    """Dual 10+6bp (16bp non-N) → [10, 6]."""
    from auto_demultiplex import get_segment_lengths
    assert get_segment_lengths('TGGGAACGCTCCGGTT') == [10, 6], \
        '16bp non-N barcode should give [10, 6]'
    print('PASSED')
    return True

def test_unit_get_segment_lengths_n_at_end():
    """First 10 non-N, last 10 N → [10, 10]."""
    from auto_demultiplex import get_segment_lengths
    assert get_segment_lengths('TGGGAACGCTNNNNNNNNNN') == [10, 10], \
        'First-10-real last-10-N barcode should give [10, 10]'
    print('PASSED')
    return True

def test_unit_get_segment_lengths_n_at_start():
    """First 10 N, last 10 non-N → [10, 10]."""
    from auto_demultiplex import get_segment_lengths
    assert get_segment_lengths('NNNNNNNNNNCCGAGACGTT') == [10, 10], \
        'First-10-N last-10-real barcode should give [10, 10]'
    print('PASSED')
    return True

def test_unit_get_segment_lengths_6bp_barcode():
    """6bp barcode followed by 14 N's → [6, 14]."""
    from auto_demultiplex import get_segment_lengths
    assert get_segment_lengths('CCGGATNNNNNNNNNNNNNN') == [6, 14], \
        '6bp real + 14 N barcode should give [6, 14]'
    print('PASSED')
    return True

def test_unit_segment_rc_single():
    """Single-segment RC: RC of whole extracted sequence."""
    from auto_demultiplex import apply_segment_rc
    # AGCGTTCCCA is RC of TGGGAACGCT
    assert apply_segment_rc('AGCGTTCCCA', [10]) == 'TGGGAACGCT'
    print('PASSED')
    return True

def test_unit_segment_rc_dual():
    """Dual-segment RC: each half RC'd independently, then concatenated."""
    from auto_demultiplex import apply_segment_rc
    # R2 last 20bp = AGCGTTCCCA + AACGTCTCGG
    # RC each: TGGGAACGCT + CCGAGACGTT
    result = apply_segment_rc('AGCGTTCCCAAACGTCTCGG', [10, 10])
    assert result == 'TGGGAACGCTCCGAGACGTT', f'got {result}'
    print('PASSED')
    return True

def test_unit_segment_rc_6bp():
    """6bp barcode segment RC: only first 6bp RC'd."""
    from auto_demultiplex import apply_segment_rc
    # ATCCGG is RC of CCGGAT
    result = apply_segment_rc('ATCCGGTAGCAAAAAAAAAA', [6, 14])
    assert result[:6] == 'CCGGAT', f'first 6 should be CCGGAT, got {result[:6]}'
    print('PASSED')
    return True


# ─── INTEGRATION TESTS ───────────────────────────────────────────────────────

def test_case1_auto():
    print('\n=== Test Case 1: PE Single Barcode Auto-detect ===')
    output_dir = '/tmp/test_case1'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case1_10r_1.fq.gz',
           '-2', 'examples/Cases/case1_10r_2.fq.gz',
           '-o', output_dir]

    if not run_cmd(cmd):
        return False

    total = count_reads(os.path.join(output_dir, 'barcode1_1.fq.gz')) * 5
    if total != 25000:
        print(f'ERROR: Expected 25000 reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case1_provided():
    print('\n=== Test Case 1: Provided Barcode ===')
    output_dir = '/tmp/test_case1_provided'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case1_10r_1.fq.gz',
           '-2', 'examples/Cases/case1_10r_2.fq.gz',
           '-B', 'examples/Cases/case1.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    total = count_reads(os.path.join(output_dir, 'barcode1_1.fq.gz')) * 5
    if total != 25000:
        print(f'ERROR: Expected 25000 reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case2_provided():
    print('\n=== Test Case 2: Dual 10+10 Barcode (provided) ===')
    output_dir = '/tmp/test_case2_provided'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case2_10r_10r_1.fq.gz',
           '-2', 'examples/Cases/case2_10r_10r_2.fq.gz',
           '-B', 'examples/Cases/case2.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    barcodes = ['barcode1', 'barcode2', 'barcode3', 'barcode4', 'barcode5']
    total = total_demuxed(output_dir, barcodes)
    if total != 25000:
        print(f'ERROR: Expected 25000 demuxed reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case3_provided():
    print('\n=== Test Case 3: Dual 10+6 Barcode (provided) ===')
    output_dir = '/tmp/test_case3_provided'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case3_10r_6r_1.fq.gz',
           '-2', 'examples/Cases/case3_10r_6r_2.fq.gz',
           '-B', 'examples/Cases/case3.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    barcodes = ['barcode1', 'barcode2', 'barcode3', 'barcode4', 'barcode5']
    total = total_demuxed(output_dir, barcodes)
    if total != 25000:
        print(f'ERROR: Expected 25000 demuxed reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case4_n_masked():
    print('\n=== Test Case 4: First barcode only, N-masked second (provided) ===')
    output_dir = '/tmp/test_case4_provided'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case4_10r_10n_1.fq.gz',
           '-2', 'examples/Cases/case4_10r_10n_2.fq.gz',
           '-B', 'examples/Cases/case4.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    barcodes = ['barcode1', 'barcode2', 'barcode3', 'barcode4', 'barcode5']
    total = total_demuxed(output_dir, barcodes)
    if total != 25000:
        print(f'ERROR: Expected 25000 demuxed reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case5_n_masked():
    print('\n=== Test Case 5: Second barcode only, N-masked first (provided) ===')
    output_dir = '/tmp/test_case5_provided'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case5_10n_10r_1.fq.gz',
           '-2', 'examples/Cases/case5_10n_10r_2.fq.gz',
           '-B', 'examples/Cases/case5.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    barcodes = ['barcode1', 'barcode2', 'barcode3', 'barcode4', 'barcode5']
    total = total_demuxed(output_dir, barcodes)
    if total != 25000:
        print(f'ERROR: Expected 25000 demuxed reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case6_provided():
    print('\n=== Test Case 6: Dual barcode (provided) ===')
    output_dir = '/tmp/test_case6_provided'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case6_10r_10f_1.fq.gz',
           '-2', 'examples/Cases/case6_10r_10f_2.fq.gz',
           '-B', 'examples/Cases/case6.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    barcodes = ['barcode1', 'barcode2', 'barcode3', 'barcode4', 'barcode5']
    total = total_demuxed(output_dir, barcodes)
    if total != 25000:
        print(f'ERROR: Expected 25000 demuxed reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case7_mixed():
    print('\n=== Test Case 7: Mixed single+dual (provided) ===')
    output_dir = '/tmp/test_case7_provided'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case7_10n_10r_singleAndDual_1.fq.gz',
           '-2', 'examples/Cases/case7_10n_10r_singleAndDual_2.fq.gz',
           '-B', 'examples/Cases/case7.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    barcodes = [f'barcode{i}' for i in range(1, 11)]
    total = total_demuxed(output_dir, barcodes)
    if total != 100000:
        print(f'ERROR: Expected 100000 demuxed reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case8_dual():
    print('\n=== Test Case 8: Dual 10+10 barcodes (provided) ===')
    output_dir = '/tmp/test_case8_dual'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case8_dual10_single6_1.fq.gz',
           '-2', 'examples/Cases/case8_dual10_single6_2.fq.gz',
           '-B', 'examples/Cases/case8.dual.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    barcodes = ['barcode1', 'barcode2', 'barcode3', 'barcode4', 'barcode5']
    total = total_demuxed(output_dir, barcodes)
    if total != 50000:
        print(f'ERROR: Expected 50000 demuxed reads, got {total}')
        return False

    print('PASSED')
    return True

def test_case8_single():
    print('\n=== Test Case 8: Single 6bp barcodes (provided) ===')
    output_dir = '/tmp/test_case8_single'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    # MGI recommends 0 mismatches for 6bp barcodes; without it, dual-barcode
    # reads in the same file can falsely match 6bp single barcodes at distance 2.
    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case8_dual10_single6_1.fq.gz',
           '-2', 'examples/Cases/case8_dual10_single6_2.fq.gz',
           '-B', 'examples/Cases/case8.single.barcode',
           '-r', '-m', '0', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    barcodes = [f'barcode{i}' for i in range(6, 11)]
    total = total_demuxed(output_dir, barcodes)
    if total != 50000:
        print(f'ERROR: Expected 50000 demuxed reads, got {total}')
        return False

    print('PASSED')
    return True

def test_sequence_stat_nonzero():
    """SequenceStat.txt counts should reflect actual demuxed read counts."""
    print('\n=== Test SequenceStat counts are non-zero ===')
    output_dir = '/tmp/test_seqstat'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case1_10r_1.fq.gz',
           '-2', 'examples/Cases/case1_10r_2.fq.gz',
           '-B', 'examples/Cases/case1.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    seq_stat = read_sequence_stat(output_dir)
    for name in ['barcode1', 'barcode2', 'barcode3', 'barcode4', 'barcode5']:
        count = seq_stat.get(name, 0)
        if count == 0:
            print(f'ERROR: SequenceStat count for {name} is 0, expected ~5000')
            return False

    print('PASSED')
    return True

def test_default_mismatch_is_2():
    """Default mismatch should be 2, matching MGI recommendation for 10bp barcodes."""
    import argparse
    # Parse auto_demultiplex.py source to verify default
    with open(os.path.join(SCRIPT_DIR, 'auto_demultiplex.py')) as f:
        source = f.read()
    if "'--mismatch'" in source or '"--mismatch"' in source:
        # Check that default=2 appears near '--mismatch'
        import re
        m = re.search(r"'--mismatch'.*?default=(\d+)", source, re.DOTALL)
        if not m:
            m = re.search(r'"--mismatch".*?default=(\d+)', source, re.DOTALL)
        if m:
            default_val = int(m.group(1))
            if default_val != 2:
                print(f'ERROR: Default mismatch is {default_val}, expected 2')
                return False
        else:
            print('ERROR: Could not parse default mismatch value')
            return False
    print('PASSED')
    return True

def test_barcode_stat_canonical_headers():
    """BarcodeStat.txt must use canonical splitBarcode column headers."""
    print('\n=== Test BarcodeStat canonical headers ===')
    output_dir = '/tmp/test_bcstat_headers'
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case1_10r_1.fq.gz',
           '-2', 'examples/Cases/case1_10r_2.fq.gz',
           '-B', 'examples/Cases/case1.barcode',
           '-r', '-o', output_dir]

    if not run_cmd(cmd):
        return False

    with open(os.path.join(output_dir, 'BarcodeStat.txt')) as f:
        header = f.readline().rstrip('\n')

    expected = '#Barcode\tCorrect\tCorrected\tTotal\tPercentage(%)'
    if header != expected:
        print(f'ERROR: Expected header:\n  {expected!r}\nGot:\n  {header!r}')
        return False

    print('PASSED')
    return True

def test_barcode_mixed_length_warns():
    """Loading a barcode file with mixed sequence lengths should print a warning."""
    print('\n=== Test mixed-length barcode file warning ===')
    import tempfile

    with tempfile.NamedTemporaryFile(mode='w', suffix='.barcode', delete=False) as f:
        f.write('bc1\tAGCGTTCCCA\n')   # 10bp
        f.write('bc2\tTAATCCTTATCATTCCTATT\n')  # 20bp
        mixed_path = f.name

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'examples/Cases/case1_10r_1.fq.gz',
           '-2', 'examples/Cases/case1_10r_2.fq.gz',
           '-B', mixed_path,
           '-r', '-o', '/tmp/test_mixed_warn']

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=SCRIPT_DIR)
    os.unlink(mixed_path)

    if 'WARNING' not in result.stderr:
        print(f'ERROR: Expected a WARNING: message in stderr about mixed barcode lengths')
        print(f'stderr: {result.stderr[:400]}')
        return False

    print('PASSED')
    return True

def test_error_handling():
    print('\n=== Test Error Handling ===')

    cmd = ['python3', 'auto_demultiplex.py',
           '-1', 'nonexistent.fq.gz',
           '-o', '/tmp/test_error']

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=SCRIPT_DIR)
    if result.returncode == 0:
        print('ERROR: Should have failed for nonexistent file')
        return False

    if 'not found' not in result.stderr.lower():
        print('ERROR: Should show "not found" error')
        return False

    print('PASSED')
    return True


def main():
    unit_tests = [
        ('hamming N-wildcard', test_unit_hamming_n_wildcard),
        ('hamming normal', test_unit_hamming_normal),
        ('segment_lengths single 10bp', test_unit_get_segment_lengths_single),
        ('segment_lengths dual 10+10', test_unit_get_segment_lengths_dual_10_10),
        ('segment_lengths dual 10+6', test_unit_get_segment_lengths_dual_10_6),
        ('segment_lengths N-at-end', test_unit_get_segment_lengths_n_at_end),
        ('segment_lengths N-at-start', test_unit_get_segment_lengths_n_at_start),
        ('segment_lengths 6bp barcode', test_unit_get_segment_lengths_6bp_barcode),
        ('segment_rc single', test_unit_segment_rc_single),
        ('segment_rc dual', test_unit_segment_rc_dual),
        ('segment_rc 6bp', test_unit_segment_rc_6bp),
    ]

    integration_tests = [
        ('Case 1 Auto-detect', test_case1_auto),
        ('Case 1 Provided Barcode', test_case1_provided),
        ('Case 2 Dual 10+10 Provided', test_case2_provided),
        ('Case 3 Dual 10+6 Provided', test_case3_provided),
        ('Case 4 N-masked second', test_case4_n_masked),
        ('Case 5 N-masked first', test_case5_n_masked),
        ('Case 6 Dual provided', test_case6_provided),
        ('Case 7 Mixed single+dual', test_case7_mixed),
        ('Case 8 Dual 10+10', test_case8_dual),
        ('Case 8 Single 6bp', test_case8_single),
        ('SequenceStat non-zero', test_sequence_stat_nonzero),
        ('Default mismatch=2', test_default_mismatch_is_2),
        ('BarcodeStat canonical headers', test_barcode_stat_canonical_headers),
        ('Mixed-length barcode warning', test_barcode_mixed_length_warns),
        ('Error Handling', test_error_handling),
    ]

    all_tests = unit_tests + integration_tests

    passed = 0
    failed = 0
    failures = []

    print('\n' + '='*60)
    print('UNIT TESTS')
    print('='*60)
    for name, test_fn in unit_tests:
        try:
            print(f'\n--- {name} ---')
            if test_fn():
                passed += 1
            else:
                failed += 1
                failures.append(name)
        except Exception as e:
            print(f'EXCEPTION: {e}')
            import traceback
            traceback.print_exc()
            failed += 1
            failures.append(f'{name} (exception)')

    print('\n' + '='*60)
    print('INTEGRATION TESTS')
    print('='*60)
    for name, test_fn in integration_tests:
        try:
            if test_fn():
                passed += 1
            else:
                failed += 1
                failures.append(name)
        except Exception as e:
            print(f'EXCEPTION in {name}: {e}')
            import traceback
            traceback.print_exc()
            failed += 1
            failures.append(f'{name} (exception)')

    print(f'\n{"="*60}')
    print(f'Results: {passed} passed, {failed} failed out of {passed+failed} tests')
    if failures:
        print('FAILED:')
        for f in failures:
            print(f'  - {f}')
    return 0 if failed == 0 else 1

if __name__ == '__main__':
    sys.exit(main())
