import click
import pandas as pd
from pathlib import Path
import re
from rich.console import Console

console = Console()

def load_query_data(fasta_path: Path) -> tuple[dict, dict]:
    """Load query sequences and structures from extended FASTA format."""
    query_seq = {}
    query_struct = {}
    with open(fasta_path) as f:
        lines = [line.strip() for line in f if line.strip()]
        for i in range(0, len(lines), 3):
            name = lines[i][1:]
            query_seq[name] = lines[i + 1].lower()
            query_struct[name] = lines[i + 2]
    return query_seq, query_struct

def parse_fasta_with_structure(fasta_path: str) -> dict:
    """Parse 3-line FASTA and return dict with sequence, structure, mfe, and cai."""
    records = {}
    with open(fasta_path) as f:
        lines = [line.strip() for line in f if line.strip()]
        for i in range(0, len(lines), 3):
            header = lines[i][1:]
            seq = lines[i + 1]
            struct = lines[i + 2]
            parts = header.split()
            name = parts[0]
            mfe = float(parts[1])
            cai = float(parts[2])
            records[name] = {
                "sequence": seq,
                "structure": struct,
                "mfe": round(mfe, 1),
                "cai": cai
            }
    return records

def estimate_cds_length(seq: str) -> int:
    """Estimate CDS length by counting uppercase letters in the sequence."""
    return sum(1 for c in seq if c.isupper())

def process_sequences(
        utp_data: dict, m1psi_data: dict, len_5utr: int) -> pd.DataFrame:
    """Extract metadata and structure info from parsed FASTA dictionaries."""
    query_seq, query_struct = load_query_data(Path("ParalinearDesign/pld-tools/resources/motif_data.fasta"))
    shared_keys = utp_data.keys() & m1psi_data.keys() if m1psi_data else utp_data.keys()
    results = []
    excluded_groups = set()

    for name in sorted(shared_keys):
        utp_entry = utp_data.get(name)
        m1psi_entry = m1psi_data.get(name) if m1psi_data else {}
        if utp_entry is None:
            continue

        utp_seq = utp_entry["sequence"]
        len_cds = estimate_cds_length(utp_seq)

        if len_cds == 0:
            group = "_".join(name.split("_")[:-1])
            console.print(f"[red]CDS not detected in: {name}. Skipping group: {group}")
            excluded_groups.add(group)
            continue

        utp_struct = utp_entry["structure"]
        m1psi_struct = m1psi_entry.get("structure", "")
        len_3utr = len(utp_seq) - (len_5utr + len_cds)

        utr_motif = name.split("_")[-4]
        matched_key = next((k for k in query_seq if k in utr_motif), None)
        if matched_key is None:
            continue

        parts = name.split()[0].split("_")
        group = "_".join(parts[:-1])
        utr1 = parts[-4]
        utr2 = parts[-3]
        lam = int(parts[-2][3:]) if parts[-2].startswith("lam") else None
        typ = parts[-1]

        row = {
            "name": name,
            "group": group,
            "lambda": lam,
            "5utr": "HBB",
            "3utr-fwd": utr1,
            "3utr-bwd": utr2,
            "type": typ,
            "len-total": len(utp_seq),
            "len-cds": len_cds,
            "len-3utr": len_3utr,
            "cai": utp_entry["cai"],
            "mfe-utp": utp_entry["mfe"],
            "sequence": utp_seq,
            "structure-mfe-utp": utp_struct,
        }

        if m1psi_entry:
            row.update({
                "mfe-m1psi": m1psi_entry["mfe"],
                "structure-mfe-m1psi": m1psi_struct,
            })
        results.append(row)

    df = pd.DataFrame(results)
    if not df.empty:
        df = df[~df["group"].isin(excluded_groups)]
    return df

@click.command()
@click.option('--input-utp', required=True, type=click.Path(exists=True),
              help='UTP input FASTA file (3-line format)')
@click.option('--input-m1psi', required=False, default=None, type=click.Path(),
              help='m1ψ input FASTA file (3-line format)')
@click.option('--len-5utr', default=53, type=int, help='Length of 5\'UTR')
@click.option('--output-dir', type=str, default=".", help='Path of output directory.')
@click.option('--output-prefix', type=str, default=None, help='Prefix for output filenames.')
def main(input_utp, input_m1psi, len_5utr, output_dir, output_prefix):
    utp_data = parse_fasta_with_structure(input_utp)
    m1psi_data = parse_fasta_with_structure(input_m1psi) if input_m1psi else {}

    df = process_sequences(utp_data, m1psi_data, len_5utr)

    type_order = {"CTRL": 0, "CFLD": 1, "MASK": 2}
    df["type_order"] = df["type"].map(type_order)
    df = df.sort_values(by=["group", "lambda", "type_order"]).drop(columns="type_order")

    end_cols = [
        c for c in df.columns
        if c.startswith("sequence") or c.startswith("structure")
    ]
    new_cols = [c for c in df.columns if c not in end_cols] + end_cols
    df = df[new_cols]
    output_prefix = output_prefix or input_utp.split("/")[-1].replace(".fasta", "")
    output_path = f"{output_dir}/{output_prefix}.tsv"
    df.to_csv(output_path, sep='\t', index=False)
    print(f"=> Saved: {output_path}")

if __name__ == '__main__':
    main()
