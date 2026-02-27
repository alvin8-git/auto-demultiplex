# MGI Demultiplexing Toolkit

A self-contained toolkit for demultiplexing MGI sequencing data. Includes the official `splitBarcode` binary and an auto-detecting Python demultiplexing script.

## Features

- **splitBarcode**: Official MGI demultiplexing binary (Linux)
- **auto_demultiplex.py**: Python script with auto-detection capabilities
- **Barcode libraries**: BarcodeV2.1.txt (10bp) and BarcodeV3.0.txt (20bp)
- **Test cases**: Ready-to-use example data
- **Mismatch tolerance**: Configurable barcode matching with mismatches
- **Progress indicator**: Shows progress for large files
- **Input validation**: Checks file existence before processing

## Quick Start

### Using auto_demultiplex.py (Recommended)

```bash
# Paired-end mode (auto-detect)
python3 auto_demultiplex.py -1 examples/Cases/case1_10r_1.fq.gz \
    -2 examples/Cases/case1_10r_2.fq.gz -o output/

# Single-end mode
python3 auto_demultiplex.py -1 input.fq.gz -o output/

# With provided barcode file
python3 auto_demultiplex.py -1 r1.fq.gz -2 r2.fq.gz \
    -B BarcodeV2.1.txt -o output/
```

### Using splitBarcode (Original Binary)

```bash
# Set library path
export LD_LIBRARY_PATH=$PWD/lib:$LD_LIBRARY_PATH

# Paired-end mode
./bin/splitBarcode -B BarcodeV2.1.txt \
    -1 examples/Cases/case1_10r_1.fq.gz \
    -2 examples/Cases/case1_10r_2.fq.gz \
    -o output/ -b 200 10 1 -r

# Single-end mode
./bin/splitBarcode -B BarcodeV2.1.txt \
    -1 examples/Cases/case1.barcode \
    -o output/ -b 100 10 1
```

## Installation

### Prerequisites

- Linux (for splitBarcode binary)
- Python 3.6+

### Setup

1. Clone or download this repository
2. Navigate to the directory:
   ```bash
   cd auto_demultiplex
   ```

That's it! Everything is self-contained.

## Usage

### auto_demultiplex.py Options

| Option | Description |
|--------|-------------|
| `-1 FILE` | Input FASTQ file (R1 for PE) |
| `-2 FILE` | Input FASTQ file R2 (for paired-end) |
| `-B FILE` | Barcode file (tsv: name\\tsequence) |
| `-o DIR` | Output directory |
| `-r` | Apply reverse complement |
| `--sample_size N` | Reads to sample for auto-detection (default: 50000) |

### splitBarcode Options

| Option | Description |
|--------|-------------|
| `-B FILE` | Barcode file |
| `-1 FILE` | Input FASTQ |
| `-2 FILE` | R2 FASTQ for PE mode |
| `-o DIR` | Output directory |
| `-b READLEN BCLEN MIS` | Position: read_len, barcode_len, mismatches |
| `-r` | Reverse complement barcodes |

## Test Cases

### Case 1: Single 10bp Barcode (PE)
```bash
python3 auto_demultiplex.py -1 examples/Cases/case1_10r_1.fq.gz \
    -2 examples/Cases/case1_10r_2.fq.gz -o test_output/
```

### Case 2: Dual 10bp+10bp Barcode (PE)
```bash
python3 auto_demultiplex.py -1 examples/Cases/case2_10r_10r_1.fq.gz \
    -2 examples/Cases/case2_10r_10r_2.fq.gz -o test_output/
```

### Case 3: Dual 10bp+6bp Barcode (PE)
```bash
python3 auto_demultiplex.py -1 examples/Cases/case3_10r_6r_1.fq.gz \
    -2 examples/Cases/case3_10r_6r_2.fq.gz -o test_output/
```

### Verify Results

```bash
python3 checkSplitResult.py test_output/
```

### Run Test Suite

```bash
python3 test_auto_demultiplex.py
```

This runs automated tests for all test cases and verifies error handling.

## Output Files

- `{barcode_id}_1.fq.gz` - R1 reads per barcode
- `{barcode_id}_2.fq.gz` - R2 reads per barcode (PE mode)
- `ambiguous_1.fq.gz` / `ambiguous_2.fq.gz` - Unmatched reads
- `BarcodeStat.txt` - Demultiplexing statistics
- `SequenceStat.txt` - Barcode distribution
- `used_barcodes.txt` - Barcodes used

## Barcode File Format

Tab-separated file:
```
barcode1	TGGGAACGCT
barcode2	TAATCCTTAT
barcode3	CAATGCGTAG
```

## MGI Sequencing Chemistry

- **PE mode**: Barcodes at end of read 2 (R2), reverse complemented
- **SE mode**: Barcodes at end of read 1 (R1)
- Common lengths: 10bp (single), 20bp (dual)

## Directory Structure

```
auto_demultiplex/
├── bin/
│   └── splitBarcode          # Main executable
├── lib/                       # Shared libraries
├── examples/
│   └── Cases/                 # Test data
├── BarcodeV2.1.txt           # 10bp barcodes
├── BarcodeV3.0.txt           # 20bp barcodes
├── auto_demultiplex.py       # Auto-detecting script
├── checkSplitResult.py       # Validation script
└── README.md
```

## License

This toolkit includes pre-compiled binaries from MGI/BGI. Please refer to their licensing terms.

## References

- splitBarcode Manual: `splitBarcode_manual.2.0.0.pdf`
- SOP: `SOP_FBS-SOP-BIN-001 Demultiplexing.pdf`
