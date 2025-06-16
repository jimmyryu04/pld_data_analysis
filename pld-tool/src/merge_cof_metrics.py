import click
import subprocess as sp
import pandas as pd
import os
import sys
from glob import glob
from pathlib import Path
from rich.console import Console

console = Console()

def write_extended_fasta(tsv_path, output_dir, filename="temp.fasta"):
    """Write 3-line FASTA with name, sequence, and structure from TSV."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    fasta_path = output_path / filename
    df = pd.read_csv(tsv_path, sep="\t")

    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            f.write(f">{row['name']} {row['mfe-utp']} {row['cai']}\n")
            f.write(f"{row['sequence']}\n")
            f.write(f"{row['structure-mfe-utp']}\n")

    return str(fasta_path)


def run_cof(cof_repo, output_dir, cof_metrics, seq_tsv_path, output_prefix):
    """Run codon-opt-factors using Snakemake and return parsed TSV result."""
    cof_repo = Path(cof_repo)
    if not cof_repo.exists():
        console.print(f"[red]Cannot find codon-opt-factors directory at {cof_repo}")
        sys.exit(1)

    try:
        sp.run(["git", "-C", str(cof_repo), "checkout", "PLD"], check=True)
    except sp.CalledProcessError:
        console.print("[red]Failed to switch to PLD branch in codon-opt-factors")
        sys.exit(1)

    fasta_path = write_extended_fasta(seq_tsv_path, output_dir)
    console.print(f"[green]Extended FASTA saved to: {fasta_path}")

    metrics_line = f"METRICS: {str(cof_metrics)}\n" if cof_metrics else ""

    config_path = Path(output_dir) / "config.yaml"
    with open(config_path, "w") as f:
        f.write(
            f"{metrics_line}EXTENDED_FASTA: True\n"
            f"EXTENDED_FASTA_FILE: {os.path.abspath(fasta_path)}\n"
            f"OUTPUT_DIR: {os.path.abspath(output_dir)}"
        )

    try:
        sp.run(
            ["snakemake", "-j", "32", "--configfile", os.path.abspath(config_path), "--latency-wait", "60"],
            cwd=str(cof_repo),
            check=True,
            stderr=sp.PIPE,
            text=True,
        )
    except sp.CalledProcessError as e:
        console.print(f"[red]Snakemake failed:\n{e.stderr}")
        sys.exit(1)

    output_path = Path(output_dir) / "cof_evals_temp.tsv"
    if not output_path.exists():
        console.print("[red]Could not find cof output tsv.")
        sys.exit(1)

    try:
        os.remove(fasta_path)
        console.print(f"[dim]Deleted temporary FASTA: {fasta_path}")
    except Exception as e:
        console.print(f"[yellow]Warning: Failed to delete {fasta_path}: {e}")
    try:
        os.remove(config_path)
        console.print(f"[dim]Deleted temporary config: {config_path}")
    except Exception as e:
        console.print(f"[yellow]Warning: Failed to delete {config_path}: {e}")

    return pd.read_csv(output_path, sep="\t")


def load_cof_tsv(cof_tsv_path, seq_tsv_path, cof_repo, output_dir, cof_metrics, output_prefix):
    """Load existing COF TSV or run COF pipeline."""
    if cof_tsv_path and os.path.exists(cof_tsv_path):
        return pd.read_csv(cof_tsv_path, sep="\t")
    else:
        console.print("[yellow]Input cof TSV not found. Running codon-opt-factors...")
        return run_cof(cof_repo, output_dir, cof_metrics, seq_tsv_path, output_prefix)


def simplify_folder_name(folder: str) -> str:
    """Simplify folder names for easier reading."""
    return (
        folder
        .replace("ViennaRNA", "V")
        .replace("partition", "p")
        .replace("fold", "f")
        .replace("m1psi", "m1")
        .replace("LinearFold", "LF")
        .replace("LinearPartition", "LP")
        .replace(":", "")
    )

def merge_cof_with_seq(cof_df, seq_df, merge_metrics=None):
    """Merge COF metrics with sequence data."""
    cof_df = cof_df[
        (cof_df["variant"] == "full") &
        (~cof_df["folder"].str.contains("LinearPartition|Linearfold", na=False))
    ].fillna("")

    if merge_metrics is not None:
        cof_df = cof_df[cof_df["metric"].isin(merge_metrics)]

    cof_df["folder"] = cof_df["folder"].apply(simplify_folder_name)
    cof_df["metric_label"] = cof_df.apply(
        lambda row: row["metric"] if row["folder"] == "" else f"{row['metric']}_{row['folder']}",
        axis=1
    )

    cof_wide = cof_df.pivot_table(
        index="seqname", columns="metric_label", values="value", aggfunc="first"
    ).reset_index()
    cof_wide.rename(columns={"seqname": "name"}, inplace=True)

    merged_df = pd.merge(seq_df, cof_wide, on="name", how="left")

    end_cols = [c for c in merged_df.columns if c.startswith("sequence") or c.startswith("structure")]
    new_cols = [c for c in merged_df.columns if c not in end_cols] + end_cols
    return merged_df[new_cols]


@click.command()
@click.option("--input-utp-tsv", "seq_tsv", required=True, type=click.Path(exists=True),
              help="Path to input UTP TSV file.")
@click.option("--output-dir", default=".", show_default=True, type=click.Path(),
              help="Directory to save the merged TSV.")
@click.option("--output-prefix", default="cofm", help="Prefix for output TSV file.")
@click.option("--cof-tsv", default=None, type=click.Path(),
              help="Optional: precomputed COF TSV file.")
@click.option("--cof-repo", default=None, type=click.Path(),
              help="Path to the codon-opt-factors repo (required if running COF).")
@click.option("--cof-metrics", multiple=True,
              help="Metrics to evaluate with codon-opt-factors.")
@click.option("--merge-metrics", multiple=True,
              help="Subset of COF metrics to keep in merged output.")
def main(seq_tsv, output_dir, output_prefix, cof_tsv, cof_repo, cof_metrics, merge_metrics):
    """Merge COF metrics into UTP sequence TSV."""
    output_tsv_path = os.path.join(output_dir, output_prefix + "_merged.tsv")
    os.makedirs(output_dir, exist_ok=True)

    seq_df = pd.read_csv(seq_tsv, sep="\t")
    cof_df = load_cof_tsv(cof_tsv, seq_tsv, cof_repo, output_dir, cof_metrics or None, output_prefix)
    merged_df = merge_cof_with_seq(cof_df, seq_df, merge_metrics or None)

    merged_df.to_csv(output_tsv_path, sep="\t", index=False)
    console.print(f"[bold green]=> Merged TSV saved to: {output_tsv_path}")

    if cof_repo:
        eval_tsv = Path(output_dir) / "cof_evals_temp.tsv"
        newname = Path(output_dir) / f"{output_prefix}_cof_evals.tsv"
        eval_tsv.rename(newname)
        m1psi_fasta  = Path(output_dir) / "m1psi_structure_temp.fasta"
        newname = Path(output_dir) / f"{output_prefix}_m1psi_structure.fasta"
        m1psi_fasta.rename(newname)

if __name__ == "__main__":
    main()
