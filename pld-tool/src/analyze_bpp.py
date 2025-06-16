import RNA
import numpy as np
import pandas as pd
import click
from pathlib import Path
from rich.console import Console
from rich.progress import (
    Progress, BarColumn, TextColumn,
    TimeElapsedColumn, TimeRemainingColumn,
    MofNCompleteColumn
)
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Optional, Tuple

console = Console()

def load_query_data(fasta_path: Path) -> tuple[dict, dict]:
    """Load query sequences and structures from extended FASTA format."""
    query_seq = {}
    query_struct = {}

    with open(fasta_path) as f:
        lines = [line.strip() for line in f if line.strip()]
        for i in range(0, len(lines), 3):
            name = lines[i][1:]
            seq = lines[i + 1].lower()
            struct = lines[i + 2]
            query_seq[name] = seq
            query_struct[name] = struct
    return query_seq, query_struct

def filter_valid_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Return only rows from groups where all sequences have the same total length."""
    valid = df.groupby("group")["len-total"].nunique()
    return df[df["group"].isin(valid[valid == 1].index)].copy()

def load_bpp_sparse(path: Path, seq_len: int) -> np.ndarray:
    """Load a sparse BPP file and return a full symmetric matrix of base pairing probabilities."""
    bpp = np.zeros((seq_len, seq_len))

    with open(path) as f:
        for line in f:
            i, j, p = line.strip().split()
            i, j = int(i) - 1, int(j) - 1
            p = float(p)
            bpp[i, j] = bpp[j, i] = p
    return bpp

def calculate_bpp_sparse(
        sequence: str, path: Path, threshold: float = 1e-6, 
        sfact: float = 1.1) -> None:
    """Calculate BPP and save (i, j, prob) if prob >= threshold."""
    sequence = sequence.upper().replace("T", "U").replace("\n", "")
    md = RNA.md(sfact=sfact)
    fc = RNA.fold_compound(sequence, md)
    struct, mfe = fc.mfe()

    fc.exp_params_rescale(mfe)
    fc.pf()
    bpp_raw = np.array(fc.bpp())[1:, 1:]
    n = bpp_raw.shape[0]

    with open(path, "w") as f:
        for i in range(n):
            for j in range(i + 1, n):
                prob = bpp_raw[i, j]
                if prob >= threshold:
                    f.write(f"{i+1}\t{j+1}\t{prob:.6f}\n")

def load_or_compute_bpp(name: str, seq: str, bpp_path: Path) -> np.ndarray:
    """Load BPP from file or compute and save if not available."""
    if bpp_path.exists() and bpp_path.stat().st_size > 0:
        console.print(f"=> [white]Path exists and not empty. Loading. ({bpp_path})")
    else:
        console.print(f"=> [green]Generating BPP for {name}")
        calculate_bpp_sparse(seq, bpp_path)
    return load_bpp_sparse(bpp_path, seq_len=len(seq))

def mean_prob_structural_elements(
    prob_matrix: np.ndarray,
    dot_bracket: str,
    region: Optional[Tuple[int, int]] = None
) -> float:
    n = len(dot_bracket)
    if region is not None:
        start, end = region
    else:
        start, end = 0, n

    if n == 0 or start == end:
        return 1.0

    pair_list = []
    unpaired_list = []
    stack = []
    for i, char in enumerate(dot_bracket):
        idx = start + i
        if char == '(':
            stack.append(idx)
        elif char == ')':
            j = stack.pop()
            if start <= i < end or start <= j < end:
                pair_list.append((j, idx))
        elif char == '.' and start <= idx < end:
            unpaired_list.append(idx)

    p_unpaired = 1.0 - prob_matrix.sum(axis=0) - prob_matrix.sum(axis=1)

    sum_of_pair_probs = 0.0
    if pair_list:
        rows, cols = zip(*pair_list)
        sum_of_pair_probs = prob_matrix[rows, cols].sum()

    unpair_probs = []
    for i in unpaired_list:
        total_pair_prob = np.sum(prob_matrix[i]) - prob_matrix[i, i]
        unpair_prob = 1.0 - total_pair_prob
        unpair_prob = np.clip(unpair_prob, 0.0, 1.0)
        unpair_probs.append(unpair_prob)

    sum_of_unpaired_probs = np.sum(unpair_probs)

    total_elements = len(pair_list) + len(unpaired_list)
    if total_elements == 0:
        return 1.0

    mean_prob = (sum_of_pair_probs + sum_of_unpaired_probs) / total_elements
    return round(mean_prob, 6)

def compute_region_bpp_avg(bpp: np.ndarray, region1: tuple, region2: tuple) -> float:
    """Return the mean BPP between two regions specified as (start, end) index tuples."""
    sub = bpp[region1[0]:region1[1], region2[0]:region2[1]]
    return round(np.mean(sub), 6) if sub.size > 0 else -1.0

def compute_bpp_metrics(row, bpp, motif_data, len_5utr):
    """Compute BPP-based metrics for a given row and its BPP matrix."""
    seq = row["sequence"].upper()
    len_cds = row["len-cds"]
    cds_start = len_5utr
    cds_end = cds_start + len_cds
    utr3_start = cds_end

    motif, struct = motif_data
    metrics = {
        "bpp-utp-cc": compute_region_bpp_avg(bpp, (cds_start, cds_end), (cds_start, cds_end)),
        "bpp-utp-c3": compute_region_bpp_avg(bpp, (cds_start, cds_end), (utr3_start, len(seq))),
        "bpp-utp-33": compute_region_bpp_avg(bpp, (utr3_start, len(seq)), (utr3_start, len(seq))),
        "bpp-utp-mm": -1.0,
        "bpp-utp-me": -1.0,
    }

    if motif and motif in seq:
        start = seq.find(motif.upper())
        end = start + len(motif)
        region = (start, end)
        metrics["bpp-utp-mm"] = round(np.mean(bpp[start:end, start:end]), 6)
        metrics["MPSE-utp"] = mean_prob_structural_elements(bpp, struct, region)
        
        nonmotif_idx = np.concatenate([np.arange(0, start), np.arange(end, len(seq))])
        metrics["bpp-utp-me"] = round(np.mean(bpp[np.ix_(np.arange(start, end), nonmotif_idx)]), 6) if len(nonmotif_idx) else -1.0
    else:
        metrics["MPSE-utp"] = -1.0
    return metrics

def compute_metrics_for_row(args):
    """Multiprocessing-safe function for computing BPP metrics from a DataFrame row."""
    row, bpp_dir, motif_dict, struct_dict, len_5utr = args
    name = row["name"]
    seq = row["sequence"].upper()
    utr = name.split("_")[1]
    key = next((k for k in motif_dict if k in utr), None)
    motif = motif_dict.get(key, None).upper() if key else None
    struct = struct_dict.get(key, None)

    bpp_path = Path(bpp_dir) / f"{name}_bpp.txt"
    bpp = load_or_compute_bpp(name, seq, bpp_path)
    return name, compute_bpp_metrics(row, bpp, (motif, struct), len_5utr)

@click.command()
@click.option("--input-tsv", required=True, type=click.Path(exists=True))
@click.option("--bpp-dir", required=True, type=click.Path())
@click.option('--output-dir', '-o', default=".", help='Output directory.')
@click.option("--len-5utr", default=53, type=int)
@click.option("--process", "-p", default=4, type=int, show_default=True, help="Number of parallel processes.")
def main(input_tsv, bpp_dir, output_dir, len_5utr, process):
    df = pd.read_csv(input_tsv, sep="\t")
    console.print(f"[cyan]Loaded input TSV: {input_tsv}")
    console.print(f"[cyan]Total rows before filtering: {len(df)}")
    df = filter_valid_groups(df)
    if df.empty:
        console.print("[red]No valid rows after filtering — nothing to process. Exiting.")
        return
    console.print(f"[bold green]Filtered to {len(df)} rows.")

    Path(bpp_dir).mkdir(parents=True, exist_ok=True)
    console.print(f"[cyan]BPP directory set to: {bpp_dir}")
    query_seq, query_struct = load_query_data(Path("/home/ryu/project/ParalinearDesign/pld-tools/resources/motif_data.fasta"))

    args_list = [
        (row, bpp_dir, query_seq, query_struct, len_5utr)
        for _, row in df.iterrows()
    ]

    results_map = {}
    with Progress(TextColumn("{task.description}"), BarColumn(), MofNCompleteColumn(),
                  TextColumn("|"), TimeElapsedColumn(), TextColumn("|"), TimeRemainingColumn(),
                  console=console) as progress:
        task = progress.add_task("Calculating and Analyzing BPP...", total=len(df))

        with ProcessPoolExecutor(max_workers=process) as executor:
            futures = {executor.submit(compute_metrics_for_row, arg): arg[0]["name"] for arg in args_list}
            for future in as_completed(futures):
                name, metrics = future.result()
                results_map[name] = metrics
                progress.update(task, advance=1)

    for col in next(iter(results_map.values())).keys():
        df[col] = [results_map[row["name"]][col] for _, row in df.iterrows()]

    end_cols = [
        c for c in df.columns
        if c.startswith("sequence") or c.startswith("structure")
    ]
    new_cols = [c for c in df.columns if c not in end_cols] + end_cols
    df = df[new_cols]

    output = f"{output_dir}/{Path(input_tsv).stem}_bpp.tsv"
    df.to_csv(output, sep="\t", index=False)
    console.print(f"[bold green]=> Saved updated TSV: {output}")
    console.print(f"[bold green]Completed BPP metric calculation for {len(df)} sequences. ✅")

if __name__ == "__main__":
    main()
