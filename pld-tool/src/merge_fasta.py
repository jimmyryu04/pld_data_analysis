import re
import click
from pathlib import Path
from typing import List, Dict, Tuple

# 1-1. 첫 번째 UTR 목록 및 순서 (K 시리즈)
FWD_UTR_LIST = [
    'K4', 'eK4', 'mK4',
    'K5', 'eK5', 'mK5',
    'K3', 'mK3'
]
FWD_UTR_ORDER_MAP = {utr: i for i, utr in enumerate(FWD_UTR_LIST)}

# 1-2. 두 번째 UTR 목록 및 순서 (HBA/HBB 시리즈)
BWD_UTR_LIST = [
    'only', 'HBA30', 'HBA45', 'HBA60', 'HBA75',
    'HBB30', 'HBB45', 'HBB60', 'HBB75'
]
BWD_UTR_ORDER_MAP = {utr: i for i, utr in enumerate(BWD_UTR_LIST)}


def parse_fasta_file(file_path: Path) -> List[Dict[str, str]]:
    """
    Parses a 3-line FASTA file into a list of record dictionaries.
    """
    records = []
    if not file_path.exists():
        print(f"Error: File not found at {file_path}. Exiting.")
        exit(1)
        
    print(f"Reading records from {file_path.name}...")
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    for i in range(0, len(lines), 3):
        if i + 2 < len(lines):
            records.append({
                'header': lines[i].strip(),
                'sequence': lines[i+1].strip(),
                'structure': lines[i+2].strip()
            })
    print(f"Found {len(records)} records.")
    return records

def get_sort_key(record: Dict[str, str]) -> Tuple[int, int, int, int, int]:
    """
    Creates a multi-level sort key based on the final user-defined hierarchy.
    Sorts by: 1. BWD UTR, 2. FWD UTR, 3. Lambda, 4. Sample No., 5. Type
    """
    header = record['header'].split('\t')[0]
    match = re.search(r'NLucP_([A-Za-z0-9]+)_([A-Za-z0-9]+)_lam(\d+)_(CTRL|MASK)_sample(\d+)', header)

    if not match:
        return (float('inf'), float('inf'), float('inf'), float('inf'), float('inf'))

    utr_fwd, utr_bwd, lambda_str, type_str, sample_str = match.groups()

    # 1. BWD UTR 순서
    bwd_utr_order = BWD_UTR_ORDER_MAP.get(utr_bwd, float('inf'))
    
    # 2. FWD UTR 순서
    fwd_utr_order = FWD_UTR_ORDER_MAP.get(utr_fwd, float('inf'))

    # 3. 람다 값
    lambda_num = int(lambda_str)
    
    # 4. 샘플 번호 (타입보다 먼저 정렬)
    sample_num = int(sample_str)
    
    # 5. 타입 순서 (CTRL=0, MASK=1)
    type_order = 0 if type_str == 'CTRL' else 1
    
    # 최종 정렬 키 반환 (샘플 번호와 타입 순서 변경)
    return (bwd_utr_order, fwd_utr_order, lambda_num, sample_num, type_order)


@click.command()
@click.option('--cds-file', required=True, type=click.Path(exists=True, path_type=Path), help='Path to the FASTA file containing CTRL records.')
@click.option('--mask-file', required=True, type=click.Path(exists=True, path_type=Path), help='Path to the FASTA file containing MASK records.')
@click.option('--output-file', required=True, type=click.Path(path_type=Path), help='Path for the merged, filtered, and sorted output FASTA file.')
def main(cds_file: Path, mask_file: Path, output_file: Path):
    """
    Filters FASTA files for specific UTR combinations, then merges and sorts them
    to interleave CTRL and MASK records for each sample.
    """
    ctrl_records = parse_fasta_file(cds_file)
    mask_records = parse_fasta_file(mask_file)
    all_records = ctrl_records + mask_records
    print(f"\nTotal records combined: {len(all_records)}")
    
    print(f"Filtering for FWD UTRs: {', '.join(FWD_UTR_LIST)}")
    print(f"Filtering for BWD UTRs: {', '.join(BWD_UTR_LIST)}")
    
    def get_utr_parts(header: str) -> Tuple[str, str]:
        match = re.search(r'NLucP_([A-Za-z0-9]+)_([A-Za-z0-9]+)_', header)
        return (match.group(1), match.group(2)) if match else ("", "")

    filtered_records = [
        rec for rec in all_records 
        if get_utr_parts(rec['header'])[0] in FWD_UTR_LIST and
           get_utr_parts(rec['header'])[1] in BWD_UTR_LIST
    ]
    print(f"Records after filtering: {len(filtered_records)}")
    
    print("Sorting records to interleave CTRL/MASK pairs...")
    sorted_records = sorted(filtered_records, key=get_sort_key)
    print("Sorting complete.")
    
    print(f"Writing sorted records to {output_file}...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        for record in sorted_records:
            f.write(f"{record['header']}\n")
            f.write(f"{record['sequence']}\n")
            f.write(f"{record['structure']}\n")
            
    print(f"\nSuccessfully created final interleaved file: {output_file}")
    
    print("\n--- First 20 headers of sorted output ---")
    for i, record in enumerate(sorted_records[:20]):
        print(f"{i+1:<3} {record['header']}")
    print("-----------------------------------------")

if __name__ == "__main__":
    main()