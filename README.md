# Uniqprimer (Python 3 Port)

Uniqprimer designs specific primers for bacterial sequences. It identifies genomic regions unique to a set of "include" genomes compared to a set of "exclude" genomes and designs PCR primers for those regions.

**Note:** This repository is a modernized port of the original Uniqprimer software (originally written in Python 2) to **Python 3**.

## Features
- **Python 3 compatible**
- **Conda-ready** for easy dependency management.
- **Improved Primer Design:** Uses direct `primer3_core` integration for more flexible parameters.
- **Chunking Strategy:** Handles large unique sequences by breaking them into manageable pieces.

## Installation

### 1. Using Conda (Recommended)
```bash
git clone https://github.com/dgutierrezcastillo/Uniqprimer.git
cd Uniqprimer
conda env create -f environment.yml
conda activate uniqprimer
```

### 2. Manual Installation
Ensure you have the following installed and in your PATH:
- Python >= 3.12
- MUMmer >= 3.23
- Primer3 >= 2.6.1
- EMBOSS (for primersearch)

Then install the Python package:
```bash
pip install .
```

## Usage

Once installed, you can use the `uniqprimer` command:

```bash
uniqprimer -i target_genome.fasta -x exclude_genome.fasta -o primers.txt
```

### Advanced Options
- `-i, --include`: Include FASTA file (can be specified multiple times).
- `-x, --exclude`: Exclude FASTA file (can be specified multiple times).
- `-o, --output`: Output file for primers (default: `uPrimer.txt`).
- `--chunksize`: Maximum sequence length to send to primer3 at once (default: `10000`).
- `--productsizerange`: PCR product size range (default: `200-250`).
- `--mintm`, `--maxtm`, `--opttm`: Melting temperature constraints.
- `--mingc`, `--maxgc`: GC content constraints.
- `--crossvalidate`: Verify primers against genomes using `primersearch`.

Example with relaxed parameters:
```bash
uniqprimer -i data/Xoo_KACC_10331.fasta -x data/XCCgenome.fasta --productsizerange 50-1000 --mintm 45 --maxtm 75 -o output.txt
```

## Credits
Originally developed at [Colorado State University](https://www.colostate.edu/) by the [Jan Leach Lab](http://leachlab.agsci.colostate.edu/).
Updated and modernized by Diego Gutierrez.

## License
GNU GPL v3
