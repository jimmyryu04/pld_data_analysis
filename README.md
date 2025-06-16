# ParalinearDesign: Advanced mRNA Sequence Optimization

ParalinearDesign is a command-line tool developed in Python that enhances mRNA coding sequences for improved stability and translation efficiency. It works by simultaneously optimizing two critical factors: codon usage bias (quantified by the Codon Adaptation Index) and RNA secondary structure stability (measured by Minimum Free Energy).

Building on the original [LinearDesign framework](https://doi.org/10.1038/s41586-023-06127-z), this advanced version introduces several improvements:

-   Co-folding of coding sequences (CDS) and 3' UTR regions
-   Broader sequence space for multiple candidates using diverse sampling strategies
-   Fast processing through parallel computation
-   Option to exclude specific restriction sites
-   Masking capability to prevent base-pairing between selected regions
-   Fine-tuning for exclusion of rare codons

## Core Algorithm: LinearDesign

ParalinearDesign implements and extends the core algorithm developed by Baidu Research:

> He Zhang†, Liang Zhang†, Ang Lin†, Congcong Xu†, Ziyu Li, Kaibo Liu, Boxiang Liu, Xiaopin Ma, Fanfan Zhao, Huiling Jiang, Chunxiu Chen, Haifa Shen, Hangwen Li*, David H. Mathews*, Yujian Zhang*, Liang Huang†*<sup>#</sup>. **Algorithm for Optimized mRNA Design Improves Stability and Immunogenicity.** *Nature* (2023). [DOI: 10.1038/s41586-023-06127-z](https://doi.org/10.1038/s41586-023-06127-z)
>
> († contributed equally, \* corresponding authors, <sup>#</sup> lead corresponding author)

### Licensing and Patent
-   **License:** As a derivative work, ParalinearDesign is subject to the same **non-commercial use only** license as the original LinearDesign. Please refer to `LICENSE.txt` for the full license text.
-   **Patent:** The core LinearDesign algorithm is covered by a patent filed by Baidu Research (Inventors: He Zhang, Liang Zhang, Ziyu Li, Kaibo Liu, Boxiang Liu, Liang Huang).

## Installation

### Prerequisites
-   Python 3 (version 3.9 or newer recommended)
-   A C++ compiler supporting the C++11 standard (e.g., GCC/g++ or Clang)
-   `pip` (Python package installer)

### Steps
1.  Ensure all prerequisites are installed and accessible in your system's PATH.
2.  Clone the repository:
    ```bash
    git clone https://github.com/ChangLabSNU/ParalinearDesign.git
    cd ParalinearDesign
    ```
3.  Install the package using pip. This command compiles the necessary C++ backend, installs Python dependencies, and makes the `pldesign` command available:
    ```bash
    pip install .
    ```
    *(Python dependencies installed automatically: `numpy>=2.0`, `pandas>=2.0`, `tqdm>=4.0`, `biopython>=1.5`)*

## Usage

The primary way to use ParalinearDesign is via the `pldesign` command-line interface:
```bash
pldesign --protein <AMINO_ACID_SEQUENCE> [OPTIONS]
```
Provide the target protein sequence and specify desired options.

## Command-Line Options

### Input & Output Settings
-   `-p, --protein PROTEIN` (**Required**): Target protein sequence using single-letter amino acid codes.
-   `-u, --utr3 UTR3`: Nucleotide sequence to append as the 3' Untranslated Region (UTR). (Default: No UTR)
-   `-o, --output OUTPUT`: Path to the output file. (Default: Standard output). Results are in CSV format unless `--fasta` or `--fasta-cds-only` is specified.
-   `--fasta`: Output results in FASTA format (includes CDS and 3' UTR).
-   `--fasta-cds-only`: Output results in FASTA format (includes only the Coding Sequence (CDS)).
-   `--allow-no-stop`: Allow input protein sequences that lack a standard stop codon (*, U, O).
-   `-q, --quiet`: Suppress progress bars and other non-essential console output.

### Core Optimization Parameters
-   `-l, --lambda LDLAMBDA`: Balance factor between MFE (structure stability) and CAI (codon usage). Higher values favor CAI. (Default: 1.0)
-   `-g, --gamma CAIGAMMA`: Exponent applied to CAI weights during optimization. (Default: 1.0)
-   `-c, --codonusagefile CODONUSAGEFILE`: Path to a custom codon usage frequency table. See bundled tables (e.g., `codon_usage_freq_table_human.csv`) for the required format. (Default: Bundled human table)
-   `--minimum-cai-weight WEIGHT`: Exclude codons with a CAI weight below this threshold. (Default: 0.0)
-   `--minimum-rscu RSCU`: Exclude codons with a Relative Synonymous Codon Usage (RSCU) below this threshold. (Default: 0.0)
-   `-m, --mask MASK`: Define a structural constraint to prevent base pairing between specified regions. Format: `"start5,end5,start3,end3"` (0-based indices). Can be used multiple times. *Note: Incompatible with adaptive sampling.*

### Sampling Strategies
-   `-r, --num-random-samples N`: Number of independent optimization runs (samples) to generate. (Default: 1). If N > 1, runs are parallelized.
-   `--random-rejection-rate RATE`: Probability threshold for rejection sampling during random sequence initialization. (Default: 0.01)
-   `--random-seed-start SEED`: Starting seed for the random number generator. For N samples, seeds `S` to `S+N-1` are used. (Default: 0)
-   `--dont-retry-on-failure`: If an optimization run fails for a specific seed, do not attempt to retry it.

### Adaptive Sampling (*Incompatible with `--mask`*)
This strategy iteratively refines structural constraints based on intermediate results.
-   `--adaptive-sampling N`: Enable adaptive sampling and set the number of iterations. (Default: 0, disabled)
-   `--adaptive-sampling-max-gap-size SIZE`: Max allowed gap (unpaired bases) within a stem considered for masking. (Default: 2)
-   `--adaptive-sampling-min-group-size SIZE`: Minimum length (base pairs) of a stem considered for masking. (Default: 5)
-   `--adaptive-sampling-max-group-size SIZE`: Maximum length (base pairs) of a stem considered for masking. (Default: 8)
-   `--adaptive-sampling-min-mask-regions FRAC`: Minimum fraction of sequence length to target for masking per iteration. (Default: 0.1)
-   `--adaptive-sampling-max-mask-regions FRAC`: Maximum fraction of sequence length to target for masking per iteration. (Default: 0.3)
-   `--adaptive-sampling-turnover-rate RATE`: Fraction of masks replaced between iterations (controls exploration). (Default: 0.3)

### Performance Tuning
-   `-t, --threads N`: Number of CPU threads for parallel execution when `-r > 1`. (Default: 1)

## Examples

### Basic Optimization (CSV Output)
Optimize "MNDTEAI" with defaults, print CSV results to console.
```bash
pldesign -p MNDTEAI
```
```csv
# Example Output:
sequence,structure,mfe,cai
AUGAACGAUACGGAGGCGAUCUAA,......(((.((....)))))....,-1.1,0.695
```

### Optimization with 3' UTR (FASTA Output)
Optimize "MNDTEAI", add 3' UTR "GGGAAA", output full sequence in FASTA.
```bash
pldesign -p MNDTEAI -u GGGAAA --fasta
```
```fasta
# Example Output:
>opt_1 mfe=-1.1 cai=0.695
AUGAACGAUACGGAGGCGAUCUAAGGGAAA
```

### Using Yeast Codon Table & Custom Lambda
Optimize with yeast table, emphasizing CAI more (lambda=1.5).
```bash
pldesign -p MNDTEAI -l 1.5 -c codon_usage_freq_table_yeast.csv
```

### Avoiding Restriction Sites
Optimize a sequence while excluding specified restriction sites:
```bash
pldesign -p MNDTEAI --exclude GAATTC --exclude GGATCC
```
This example avoids EcoRI and BamHI restriction sites in the optimized sequence.

### Generating Multiple Samples (Parallel)
Generate 10 samples using 4 threads, seed 100, save to `results.csv`.
```bash
pldesign -p YOUR_PROTEIN_SEQUENCE -r 10 -t 4 --random-seed-start 100 -o results.csv
```

### Using Adaptive Sampling
Perform 5 adaptive sampling iterations, save to `adaptive_results.csv`.
```bash
pldesign -p YOUR_PROTEIN_SEQUENCE --adaptive-sampling 5 -o adaptive_results.csv
```

## Authors of ParalinearDesign
This enhanced implementation (ParalinearDesign) was developed at Seoul National University by Hyeshik Chang, Jisu Lim, and Hoon Ma.
