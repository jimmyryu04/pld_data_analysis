import pandas as pd
import click
from pathlib import Path
from datetime import datetime


def filter_valid_groups(df: pd.DataFrame, log_path: Path, label: str = "") -> list:
    """Filter valid groups by structure metrics and append log to text file."""
    grouped = df.groupby("group")
    valid_groups = []

    count_missing_types = 0
    count_mfe_condition = 0
    count_bpp_c3_condition = 0
    count_bpp_mm_condition = 0

    for group, gdf in grouped:
        gdf = gdf.set_index("type")
        if not all(t in gdf.index for t in ["CTRL", "CFLD", "MASK"]):
            count_missing_types += 1
            continue

        cds = gdf.loc["CTRL"]
        utr = gdf.loc["CFLD"]
        mask = gdf.loc["MASK"]

        if isinstance(cds, pd.DataFrame):
            cds = cds.iloc[0]
        if isinstance(utr, pd.DataFrame):
            utr = utr.iloc[0]
        if isinstance(mask, pd.DataFrame):
            mask = mask.iloc[0]

        if (
            cds["mfe-utp"] <= utr["mfe-utp"] 
            or cds["mfe-utp"] <= mask["mfe-utp"]
        ):
            count_mfe_condition += 1
            continue

        if (
            cds["bpp-utp-c3"] >= utr["bpp-utp-c3"] 
            or cds["bpp-utp-c3"] >= mask["bpp-utp-c3"]
        ):
            count_bpp_c3_condition += 1
            continue

        if (
            cds["bpp-utp-ms"] > 0.8 
            or utr["bpp-utp-ms"] > 0.8 
            or cds["bpp-utp-ms"] >= mask["bpp-utp-ms"] 
            or mask["bpp-utp-ms"] < 0.8
        ):
            count_bpp_mm_condition += 1
            continue

        valid_groups.append(group)

    total_groups = len(grouped)
    valid_group_count = len(valid_groups)
    excluded_group_count = total_groups - valid_group_count

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = (
        f"[{timestamp}] label={label}\n"
        f"  total={total_groups}, remain={valid_group_count}, exclude={excluded_group_count}\n"
        f"  missing={count_missing_types}, mfe={count_mfe_condition}, "
        f"bpp_c3={count_bpp_c3_condition}, bpp_mm={count_bpp_mm_condition}\n\n"
    )
    with open(log_path, "a") as f:
        f.write(log_entry)

    return valid_groups


def save_filtered_tsv(df: pd.DataFrame, output_path: Path) -> None:
    """Save filtered DataFrame to TSV."""
    df.to_csv(output_path, sep="\t", index=False)


def save_fasta(df: pd.DataFrame, fasta_path: Path) -> None:
    """Save sequence and structure columns in FASTA (3-line) format."""
    with open(fasta_path, "w") as f:
        for _, row in df.iterrows():
            name = row["name"]
            seq = row["sequence"]
            struct = row["structure-mfe-utp"]
            f.write(f">{name}\n{seq}\n{struct}\n")


@click.command()
@click.option("--input-tsv", required=True, type=click.Path(exists=True),
              help="Path to input TSV file")
@click.option("--output-dir", default=".", type=click.Path(), help="Directory to save output files")
@click.option("--output-prefix", default=None, type=str, help="Prefix used in output filenames")
@click.option("--fasta", is_flag=True, default=False,
              help="Save filtered result in FASTA format as well")
def main(input_tsv: str, output_dir: str, output_prefix: str, fasta: bool):
    """Main CLI to filter input TSV and write filtered result and log."""
    df = pd.read_csv(input_tsv, sep="\t")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = Path("../log/filter_log.txt")
    log_path.parent.mkdir(parents=True, exist_ok=True)

    label = output_prefix or Path(input_tsv).stem
    valid_groups = filter_valid_groups(df, log_path, label)
    filtered_df = df[df["group"].isin(valid_groups)]

    if output_prefix:
        filtered_path = output_dir / f"{output_prefix}_filtered.tsv"
    else:
        filtered_path = output_dir / input_tsv.split("/")[-1].replace(".tsv", "_filtered.tsv")

    save_filtered_tsv(filtered_df, filtered_path)
    print(f"=> Filtered {len(valid_groups)} groups saved to {filtered_path}")

    if fasta:
        fasta_path = filtered_path.with_suffix(".fasta")
        save_fasta(filtered_df, fasta_path)
        print(f"=> FASTA saved to {fasta_path}")


if __name__ == "__main__":
    main()
