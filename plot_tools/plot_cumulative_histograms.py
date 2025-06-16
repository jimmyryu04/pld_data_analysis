#!/usr/bin/env python3
"""
Grid-based Cumulative Histogram Plotter for MPSE (CTRL vs MASK)
─────────────────────────────────────────────────────────────────
Generates a single PNG file per UTR family, containing a 2D grid of subplots.
Columns represent UTR types (e.g., K3, mK3), and rows represent lambda conditions.

Usage:
python plot_grid_histograms.py --input-tsv <path/to/NLucP_bpp.tsv> --output-dir <path/to/output_folder>
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import click
from pathlib import Path
import re
import numpy as np

# 분석할 메트릭 이름
METRIC = "MPSE-utp"

# UTR 패밀리 정의
UTR_FAMILIES = {
    'K3_family': ['K3', 'mK3'],
    'K4_family': ['K4', 'eK4', 'mK4'],
    'K5_family': ['K5', 'eK5', 'mK5']
}
UTR_TO_FAMILY_MAP = {utr: family for family, utrs in UTR_FAMILIES.items() for utr in utrs}


def load_and_prepare(df: pd.DataFrame) -> pd.DataFrame:
    """Parses names to create UTR types, families, and lambda values for plotting."""
    if 'name' not in df.columns:
        print("Error: 'name' column not found in TSV.", file=sys.stderr)
        return pd.DataFrame()

    def extract_info(name):
        match = re.search(r'NLucP_([A-Za-z0-9]+)_[A-Za-z0-9]+_lam(\d+)_(CTRL|MASK)_sample\d+', name)
        if match:
            utr_type, lam_str, type_str = match.groups()
            utr_family = UTR_TO_FAMILY_MAP.get(utr_type)
            return utr_family, utr_type, int(lam_str), type_str
        return None, None, None, None

    extracted_info = df['name'].apply(extract_info).apply(pd.Series)
    df[['utr_family', 'utr_type', 'lambda', 'type']] = extracted_info
    
    df = df.dropna(subset=['utr_family', 'utr_type', 'lambda', 'type', METRIC])
    df['lambda'] = df['lambda'].astype(int)
    
    return df

@click.command()
@click.option('--input-tsv', required=True, type=click.Path(exists=True, path_type=Path), help='Input TSV file (e.g., NLucP_bpp.tsv).')
@click.option('--output-dir', required=True, type=click.Path(path_type=Path), help='Directory to save the output PNG plots.')
def main(input_tsv: Path, output_dir: Path):
    """
    Main function to generate and save family-grouped 2D grid plots.
    """
    print(f"Reading records from {input_tsv.name}...")
    try:
        data = pd.read_csv(input_tsv, sep='\t')
        print(f"Found {len(data)} records.")
    except Exception as e:
        print(f"Error reading TSV file: {e}", file=sys.stderr)
        sys.exit(1)

    prepared_data = load_and_prepare(data)

    if prepared_data.empty:
        print("Error: No valid data to plot after preparation.", file=sys.stderr)
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)
    
    # UTR 패밀리별로 루프를 돌며 그림(Figure) 생성
    for family_name, family_df in prepared_data.groupby('utr_family'):
        
        # 그리드 차원 계산
        utr_types_in_family = sorted(family_df['utr_type'].unique())
        lambdas_in_family = sorted(family_df['lambda'].unique())
        
        ncols = len(utr_types_in_family)
        nrows = len(lambdas_in_family)
        
        if nrows == 0 or ncols == 0:
            continue

        print(f"Generating grid plot for {family_name} ({nrows} rows x {ncols} cols)...")

        fig, axes = plt.subplots(
            nrows=nrows, 
            ncols=ncols, 
            figsize=(7 * ncols, 6 * nrows),
            squeeze=False
        )
        
        fig.suptitle(f'Cumulative Distribution of MPSE for {family_name}', fontsize=22, weight='bold')

        # 위치 지정을 위한 맵 생성
        utr_type_to_col = {utr_type: i for i, utr_type in enumerate(utr_types_in_family)}
        lambda_to_row = {lam: i for i, lam in enumerate(lambdas_in_family)}

        # 각 서브플롯 채우기
        for (utr_type, lam), condition_df in family_df.groupby(['utr_type', 'lambda']):
            row_idx = lambda_to_row.get(lam)
            col_idx = utr_type_to_col.get(utr_type)
            
            if row_idx is None or col_idx is None:
                continue
                
            ax = axes[row_idx, col_idx]
            
            ctrl_data = condition_df[condition_df['type'] == 'CTRL'][METRIC]
            mask_data = condition_df[condition_df['type'] == 'MASK'][METRIC]

            sns.ecdfplot(data=ctrl_data, ax=ax, label=f'CTRL (n={len(ctrl_data)})', linewidth=2)
            sns.ecdfplot(data=mask_data, ax=ax, label=f'MASK (n={len(mask_data)})', linewidth=2, linestyle='--')

            ax.set_title(f'{utr_type} | λ={lam}', fontsize=14)
            ax.set_xlabel('MPSE Value', fontsize=10)
            ax.set_ylabel('Cumulative Probability', fontsize=10)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.grid(True, which='both', linestyle='--', linewidth=0.5)
            ax.legend(fontsize=10, loc='lower right')

        # 비어있는 서브플롯 숨기기
        for i in range(nrows):
            for j in range(ncols):
                # axes[i, j]에 그려진 내용이 없으면 축을 숨김
                if not axes[i, j].has_data():
                    axes[i, j].axis('off')

        fig.tight_layout(rect=[0, 0, 1, 0.96])
        output_path = output_dir / f"{family_name}_grid_plots.png"
        
        try:
            plt.savefig(output_path, dpi=300)
            print(f"  -> Grid plot saved to: {output_path}")
        except Exception as e:
            print(f"  -> Error saving plot: {e}", file=sys.stderr)
        
        plt.close(fig)

    print("\nAll plots generated successfully.")

if __name__ == '__main__':
    main()