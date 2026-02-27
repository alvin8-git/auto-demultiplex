# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run the full test suite (note: test file has hardcoded chdir to /workspace/auto_demultiplex — update line 138 if running elsewhere)
python3 test_auto_demultiplex.py

# Validate demultiplexing output (PE mode — checks paired-read name consistency)
python3 checkSplitResult.py <output_dir>

# Validate with read-count check against raw data
python3 checkSplitResult.py -r <raw_fq_dir> <output_dir>

# Run the official splitBarcode binary (requires lib/ on LD_LIBRARY_PATH)
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH
./bin/splitBarcode -B BarcodeV2.1.txt -1 r1.fq.gz -2 r2.fq.gz -o output/ -b 200 10 1 -r
```

## Architecture

The project is a self-contained toolkit with two demultiplexing paths:

**`auto_demultiplex.py`** — Python wrapper with auto-detection. The key pipeline:
1. **Detection phase** (`detect_barcode_config`): Samples `--sample_size` reads (default 50,000) from R2 (PE) or R1 (SE), extracts trailing 10bp/20bp k-mers, and matches against known barcodes in `BarcodeV2.1.txt`/`BarcodeV3.0.txt` using Hamming distance ≤ 1. Falls back to extracting the most-frequent unknown sequences if no known barcodes match.
2. **Demultiplexing phase** (`demultiplex_pe` / `demultiplex_se`): Streams through the full FASTQ, extracts the last `bc_len` bases of each R2 (PE) or R1 (SE) read, performs exact dict lookup then mismatch fallback, and writes per-barcode gzipped FASTQ output files.

**MGI chemistry conventions encoded in the script:**
- PE mode: barcode is at the **end of R2**, and is reverse-complemented relative to the reference barcode sequence
- SE mode: barcode is at the **end of R1** (no reverse complement by default)
- Barcode length is assumed to be 10bp unless a 20bp barcode library matches

**`bin/splitBarcode`** — Precompiled MGI/BGI Linux binary; requires shared libs in `lib/` via `LD_LIBRARY_PATH`. This is the reference implementation that `auto_demultiplex.py` mimics.

**`checkSplitResult.py`** — Validation utility. In PE mode, verifies all paired reads share the same read name prefix (up to `/`). Optionally compares total read count against a `.fqStat.txt` sidecar file from the raw FASTQ directory.

## Barcode files

- `BarcodeV2.1.txt`: 128 standard 10bp barcodes + barcode 999 (control spike-in)
- `BarcodeV3.0.txt`: 20bp dual-indexed barcodes (IDs in `N-N` format)
- Both use tab-separated format: `<name>\t<sequence>`
- Barcode sequences in the case `.barcode` files are **NOT** from `BarcodeV2.1.txt` — they are separate test-specific sequences

## Test cases

All 8 cases are PE, barcodes at end of R2 (reverse-complemented). R1 is always 100bp. Run with `-B <case.barcode> -r`:

| Case | R2 len | Barcode type | Expected | Script result |
|------|--------|-------------|----------|---------------|
| case1_10r | 110bp | Single 10bp | 100% | **PASS** |
| case2_10r_10r | 120bp | Dual 10+10bp | 100% | **FAIL** (0%) |
| case3_10r_6r | 116bp | Dual 10+6bp | 100% | **FAIL** (0%) |
| case4_10r_10n | 120bp | Dual 10+N (first only) | 100% | **FAIL** (0%) |
| case5_10n_10r | 120bp | Dual N+10 (second only) | 100% | **FAIL** (0%) |
| case6_10r_10f | 120bp | Dual 10r+10f (rev+fwd) | 100% | **FAIL** (0%) |
| case7_10n_10r_singleAndDual | 120bp | Mixed dual+single | 100% | **FAIL** (0%) |
| case8_dual10_single6 | 120bp | Mixed dual10+single6 | 100% | **FAIL** (0%) |

Cases 2–8 all fail because barcode file sequences are 16–20bp but the script hardcodes `barcode_length = 10`, causing `hamming_distance(10bp, 20bp)` = 20 on every read. See `Documentation.md` for the full correctness analysis.

## Known issues in auto_demultiplex.py

1. `barcode_length` hardcoded to 10 — dual barcodes completely unsupported (root cause of cases 2–8 failures)
2. No dual barcode position logic — needs two independent extractions at R2[101–110] and R2[111–120]
3. No N-wildcard handling in `hamming_distance` — N in reference should match any base
4. Default `--mismatch 1` — MGI recommends 2 for 10bp barcodes
5. SE detection scans start of reads (`seq[:bc_len]`) — SE barcodes are end-only
6. `SequenceStat.txt` counts always 0 for mismatch-corrected reads (looks up reference seq in counter keyed by read seq)
7. `test_auto_demultiplex.py` line 138: hardcoded `os.chdir('/workspace/auto_demultiplex')` — update to match environment
