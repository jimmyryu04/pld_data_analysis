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
        
        # --- 코드 수정 시작 ---
        # 원본 이름에서 파싱할 부분을 분리합니다.
        # 예: 'NLucP_mK3_HBA60_lam2_CTRL_sample1' -> ['NLucP', 'mK3', 'HBA60', 'lam2', 'CTRL', 'sample1']
        name_parts_list = name.split()[0].split("_")
        
        # 마지막 부분이 'sample'로 시작하면 파싱을 위해 임시로 제거합니다.
        if name_parts_list[-1].startswith("sample"):
            parts = name_parts_list[:-1]
        else:
            parts = name_parts_list
        # 이제 'parts' 변수는 ['NLucP', 'mK3', 'HBA60', 'lam2', 'CTRL'] 와 같은 형식을 갖게 됩니다.
        # --- 코드 수정 끝 ---

        # utr_motif 부분을 찾기 위해 원본 이름을 계속 사용하거나 수정된 parts를 사용할 수 있습니다.
        # 기존 로직은 이름의 특정 위치를 가정하므로, 수정된 parts를 사용하는 것이 더 안정적일 수 있습니다.
        # 다만, 기존 코드의 의도를 최대한 유지하기 위해 아래 로직은 그대로 둡니다.
        utr_motif = name.split("_")[-4] # 이 부분은 필요시 parts[-3] 등으로 조정될 수 있습니다.
        matched_key = next((k for k in query_seq if k in utr_motif), None)
        if matched_key is None:
            # 이 부분은 utr_motif가 올바르게 추출되지 않을 경우 문제가 될 수 있습니다.
            # 위에서 수정된 'parts' 변수를 사용하여 안정성을 높일 수 있습니다.
            # 예: utr_motif = parts[-3]
            pass # 일단 기존 로직을 유지

        group = "_".join(parts[:-1])
        utr1 = parts[-4]
        utr2 = parts[-3]
        lam = int(parts[-2][3:]) if parts[-2].startswith("lam") else None
        typ = parts[-1]

        row = {
            "name": name, # 최종 결과에는 '..._sample1'이 포함된 원본 이름을 저장
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
