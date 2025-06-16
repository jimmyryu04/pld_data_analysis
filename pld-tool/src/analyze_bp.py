import click
import pandas as pd
from typing import List, Tuple
from pathlib import Path
from rich.console import Console

console = Console()

def load_query_data(fasta_path: Path) -> Tuple[dict, dict]:
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

def get_pairs(struct: str) -> List[Tuple[int, int]]:
    """Return list of base-paired indices from dot-bracket structure."""
    stack = []
    pairs = []

    for i, c in enumerate(struct):
        if c == '(': stack.append(i)
        elif c == ')':
            if stack: pairs.append((stack.pop(), i))
    return pairs

def compute_structural_match_score(
        query_struct: str,
        query: str,
        sequence: str,
        full_structure: str) -> Tuple[List[Tuple[int, int]], int, float]:
    """Calculate match score between query motif and full structure."""
    query = query.upper().replace("U", "T")
    sequence = sequence.upper().replace("U", "T")
    start = sequence.find(query)

    if start == -1:
        return [], 0, 0.0

    query_pairs = get_pairs(query_struct)
    query_pairs_global = {(start + i, start + j) for i, j in query_pairs}
    full_pairs = set(get_pairs(full_structure))

    matched_pairs = [
        (i + 1, j + 1)
        for i, j in query_pairs_global
        if (i, j) in full_pairs or (j, i) in full_pairs
    ]
    dot_indices = [i for i, c in enumerate(query_struct) if c == '.']
    matched_dots = sum(
        1 for i in dot_indices if full_structure[start + i] == '.'
    )
    total = len(dot_indices) + len(query_pairs_global)
    satisfied = matched_dots + len(matched_pairs)
    return matched_pairs, matched_dots, round(satisfied / total, 3) if total > 0 else 0.0

def get_cds_3utr_pairs(
        structure: str, len_5utr: int, len_cds: int) -> Tuple[int, int, float]:
    """Count base pairs between CDS and 3'UTR."""
    pairs = get_pairs(structure)
    cds_start = len_5utr
    cds_end = len_5utr + len_cds - 1
    count_c3 = sum(
        1 for i, j in pairs
        if (cds_start <= i <= cds_end and j > cds_end) or
           (cds_start <= j <= cds_end and i > cds_end)
    )
    total_bp = len(pairs)
    ratio = count_c3 / total_bp if total_bp > 0 else 0
    return count_c3, total_bp, round(ratio, 4)

def count_bp_cc(
        structure: str, len_5utr: int, len_cds: int) -> int:
    """Count base pairs within CDS region."""
    pairs = get_pairs(structure)
    cds_start = len_5utr
    cds_end = len_5utr + len_cds - 1
    return sum(
        1 for i, j in pairs if cds_start <= i <= cds_end and cds_start <= j <= cds_end
    )

def count_bp_mm(
        structure: str,
        motif: str,
        query_struct: str,
        sequence: str) -> int:
    """Count base pairs within motif-motif region."""
    sequence = sequence.upper().replace("U", "T")
    motif = motif.upper().replace("U", "T")
    start = sequence.find(motif)
    
    if start == -1:
        return -1
    query_pairs = get_pairs(query_struct)
    global_pairs = {(start + i, start + j) for i, j in query_pairs}
    full_pairs = set(get_pairs(structure))
    return sum((i, j) in full_pairs or (j, i) in full_pairs for (i, j) in global_pairs)

def filter_valid_groups(df: pd.DataFrame) -> pd.DataFrame:
    """Filter groups with consistent sequence lengths and expected total length."""
    df["len-total"] = df["sequence"].str.len()
    valid_groups_1 = df.groupby("group")["len-total"].nunique()
    valid_groups_1 = valid_groups_1[valid_groups_1 == 1].index

    df["expected_total"] = 53 + df["len-cds"] + df["len-3utr"]
    group_match = df["len-total"] == df["expected_total"]
    valid_groups_2 = df[group_match].groupby("group").size()
    all_rows_per_group = df.groupby("group").size()

    valid_groups_2, all_rows_per_group = valid_groups_2.align(
        all_rows_per_group, join="inner"
    )
    fully_valid_groups = valid_groups_2[valid_groups_2 == all_rows_per_group].index
    valid_groups = valid_groups_1.intersection(fully_valid_groups)
    df_filtered = df[df["group"].isin(valid_groups)].copy()

    df_filtered.drop(columns=["expected_total"], inplace=True)
    console.print(f"[bold green] Filtered to {len(df_filtered) // 3} groups.")
    return df_filtered

def compute_structure_metrics(
        df: pd.DataFrame, len_5utr: int) -> pd.DataFrame:
    """Compute structural metrics for all rows in the dataframe."""
    query_seq, query_struct = load_query_data(Path("ParalinearDesign/pld-tools/resources/motif_data.fasta"))

    bp_utp_all_list = []
    bp_utp_cc_list = []
    bp_utp_c3_list = []
    bp_utp_mm_list = []
    bp_utp_ms_list = []

    for _, row in df.iterrows():
        key = row['3utr-fwd']
        sequence = row['sequence']
        len_cds = row['len-cds']
        struct_utp = row.get('structure-mfe-utp', '')
        matched_key = next((k for k in query_seq if k in key), None)

        if matched_key and struct_utp:
            bp_all = len(get_pairs(struct_utp))
            bp_cc = count_bp_cc(struct_utp, len_5utr, len_cds)
            bp_c3, _, _ = get_cds_3utr_pairs(struct_utp, len_5utr, len_cds)
            motif = query_seq[matched_key]
            q_struct = query_struct[matched_key]
            bp_mm = count_bp_mm(struct_utp, motif, q_struct, sequence)
            _, _, bp_ms = compute_structural_match_score(
                q_struct, motif, sequence, struct_utp
            )
        else:
            bp_all = bp_cc = bp_c3 = bp_mm = bp_ms = -1

        bp_utp_all_list.append(bp_all)
        bp_utp_cc_list.append(bp_cc)
        bp_utp_c3_list.append(bp_c3)
        bp_utp_mm_list.append(bp_mm)
        bp_utp_ms_list.append(bp_ms)

    df["bp-utp-all"] = bp_utp_all_list
    df["bp-utp-cc"] = bp_utp_cc_list
    df["bp-utp-c3"] = bp_utp_c3_list
    df["bp-utp-mm"] = bp_utp_mm_list
    df["bp-utp-ms"] = bp_utp_ms_list

    end_cols = [
        c for c in df.columns
        if c.startswith("sequence") or c.startswith("structure")
    ]
    new_cols = [c for c in df.columns if c not in end_cols] + end_cols
    return df[new_cols]

@click.command()
@click.option('--input-tsv', '-i', multiple=True, required=True)
@click.option('--output-dir', '-o', default=".", help='Output directory.')
@click.option('--len-5utr', default=53, help='Length of 5′UTR.')
def main(input_tsv, output_dir, len_5utr):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for tsv_path in input_tsv:
        df = pd.read_csv(tsv_path, sep='\t')
        df = filter_valid_groups(df)
        df = compute_structure_metrics(df, len_5utr)
        df.to_csv(tsv_path, sep='\t', index=False)
        print(f"Updated and saved: {tsv_path}")

if __name__ == '__main__':
    main()
