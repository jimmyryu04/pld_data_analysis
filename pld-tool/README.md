# pld-tools
This repository contains scripts for repeatedly running ParalinearDesign and analyzing its results, intended to be used as a submodule.

### run_pld.py

| Option            | Type    | Required | Default   | Description                                                       |
| ----------------- | ------- | -------- | --------- | ----------------------------------------------------------------- |
| `--utr3-fwd`      | Path    | ✅ Yes    | –         | FASTA file with forward 3′ UTR candidate sequences                |
| `--utr3-bwd`      | Path    | ✅ Yes    | –         | FASTA file with backward 3′ UTR candidate sequences               |
| `--cds-fasta`     | Path    | ✅ Yes    | –         | FASTA file with amino acid sequences for CDS generation           |
| `--output-dir`    | String  | ❌ No     | `"."`     | Directory to save output files and logs                           |
| `--output-prefix` | String  | ❌ No     | `""`      | Prefix to prepend to output filenames                             |
| `--cds-only`      | Flag    | ❌ No     | `False`   | Run CDS-only optimization (no UTR context)                        |
| `--with-utr`      | Flag    | ❌ No     | `False`   | Include 3′UTR context during CDS optimization                     |
| `--with-mask`     | Flag    | ❌ No     | `False`   | Apply motif masking during 3′UTR-aware optimization               |
| `--lambda`        | String  | ❌ No     | `"2,3,4"` | Comma-separated list of lambda values (e.g., loop weight factors) |
| `--process`, `-p` | Integer | ❌ No     | `4`       | Number of parallel threads/processes to use                       |

---
### parse_pld_result.py

| Option            | Type    | Required | Default             | Description                                                         |
| ----------------- | ------- | -------- | ------------------- | ------------------------------------------------------------------- |
| `--input-utp`     | Path    | ✅ Yes    | –                   | Input UTP FASTA file (3-line format with name, sequence, structure) |
| `--input-m1psi`   | Path    | ❌ No     | `None`              | Optional m1ψ FASTA file (3-line format)                             |
| `--len-5utr`      | Integer | ❌ No     | `53`                | Length of the 5′UTR region (used to infer CDS and 3′UTR lengths)    |
| `--output-dir`    | String  | ❌ No     | `"."`               | Output directory to save the result `.tsv` file                     |
| `--output-prefix` | String  | ❌ No     | from input filename | Prefix for output `.tsv` file (default: derived from `--input-utp`) |

---
### analyze_bp.py

| Option               | Type    | Required | Default | Description                                                              |
| -------------------- | ------- | -------- | ------- | ------------------------------------------------------------------------ |
| `--input-tsv`, `-i`  | Path(s) | ✅ Yes    | –       | One or more input TSV file(s) with sequence and structure information    |
| `--output-dir`, `-o` | Path    | ❌ No     | `"."`   | Directory to save the updated TSV files                                  |
| `--len-5utr`         | Integer | ❌ No     | `53`    | Length of the 5′UTR region (used to separate CDS and 3′UTR for analysis) |

---
### analyze_bpp.py

| Option               | Type    | Required | Default | Description                                                        |
| -------------------- | ------- | -------- | ------- | ------------------------------------------------------------------ |
| `--input-tsv`        | Path    | ✅ Yes    | –       | Input TSV file containing sequence and metadata                    |
| `--bpp-dir`          | Path    | ✅ Yes    | –       | Directory to load or save base-pair probability (BPP) `.txt` files |
| `--output-dir`, `-o` | Path    | ❌ No     | `"."`   | Directory to save the updated output TSV file                      |
| `--len-5utr`         | Integer | ❌ No     | `53`    | Length of the 5′ UTR (used to partition CDS and 3′UTR regions)     |
| `--process`, `-p`    | Integer | ❌ No     | `4`     | Number of parallel processes                                       |

---
### merge_cof_metrics.py

| Option            | Type   | Required | Default    | Description                                                                    |
| ----------------- | ------ | -------- | ---------- | ------------------------------------------------------------------------------ |
| `--input-utp-tsv` | Path   | ✅ Yes    | –          | Input TSV file containing UTP sequence and structure info                      |
| `--output-dir`    | String | ❌ No     | `"."`      | Directory to save the final merged TSV file                                    |
| `--output-prefix` | String | ❌ No     | `"cofm"`   | Prefix for output filename                                                     |
| `--cof-tsv`       | Path   | ❌ No     | `None`     | Optional precomputed COF TSV file; if not provided, `codon-opt-factors` is run |
| `--cof-repo`      | Path   | ❌ No     | `None`     | Path to the `codon-opt-factors` repo (required if running COF)                 |
| `--cof-metrics`   | List   | ❌ No     | all        | Metrics to evaluate using `codon-opt-factors` (e.g., `start_structure aup`)    |
| `--merge-metrics` | List   | ❌ No     | all        | Subset of metrics to keep when merging (e.g., `start_structure aup_u3`)        |

---
### filter.py

| Option            | Type   | Required | Default                     | Description                                                  |
| ----------------- | ------ | -------- | --------------------------- | ------------------------------------------------------------ |
| `--input-tsv`     | Path   | ✅ Yes    | –                           | Input TSV file containing sequence and structure information |
| `--output-dir`    | Path   | ❌ No     | `"."`                       | Directory to save filtered TSV and optional FASTA output     |
| `--output-prefix` | String | ❌ No     | derived from input filename | Prefix for output file naming                                |
| `--fasta`         | Flag   | ❌ No     | `False`                     | If set, save filtered result in 3-line FASTA format          |

