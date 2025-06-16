#!/usr/bin/env python3
"""
Aesthetic scatter plot: MPSE (CTRL vs MASK) with filtering and integrated statistics.
───────────────────────────────────────────────────────────────────────────────────
python scatter_plot_mpse_random_sampling_filtering.py --input-file <input.tsv> --output-file <out.png>
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import click
from pathlib import Path

# 분석할 메트릭 이름
METRIC = "MPSE-utp"
# 필터링 임계값
THRESHOLD = 0.7

def load_and_prepare(df: pd.DataFrame) -> pd.DataFrame:
    """Loads and prepares data for plotting, handling random sampling names."""
    if 'name' not in df.columns:
        print("Error: 'name' column not found in TSV.", file=sys.stderr)
        return pd.DataFrame()

    df['type'] = df['name'].str.extract(r'_(CTRL|MASK)')
    df['group'] = df['name'].str.replace(r'_(CTRL|MASK)', '', regex=True)
    
    sub = df[df["type"].isin(["CTRL", "MASK"])].dropna(subset=['type', 'group', METRIC])
    
    if sub.empty:
        return pd.DataFrame()

    sub = sub.drop_duplicates(subset=['group', 'type'])
    pvt = sub.pivot(index="group", columns="type", values=METRIC).dropna()
    
    # 이 부분은 사용자님의 원본 로직을 그대로 유지합니다.
    # UTR 정보가 필요한 경우에 대한 처리입니다.
    if "3utr-fwd" in df.columns:
        utr_map_df = df[df["type"] == "CTRL"].drop_duplicates(subset=['group'])
        utr_map = utr_map_df.set_index("group")["3utr-fwd"]
        pvt = pvt.join(utr_map)
    else:
        # 3utr-fwd 컬럼이 없는 경우를 대비한 처리
        pvt["3utr-fwd"] = "N/A"
        
    pvt["filtered"] = (pvt["CTRL"] > THRESHOLD) | (pvt["MASK"] > THRESHOLD)
    return pvt

def calculate_stats(df: pd.DataFrame) -> dict:
    """Calculates statistics based on the filtering threshold."""
    total_points = len(df)
    if total_points == 0:
        return {}
    
    low_both = df[(df["CTRL"] < THRESHOLD) & (df["MASK"] < THRESHOLD)]
    filtered_count = len(df[df["filtered"]])
    
    stats = {
        "Total points": total_points,
        "Filtered (MPSE > 0.7)": filtered_count,
        "Retained (MPSE < 0.7)": len(low_both),
        "Percentage filtered": f"{(filtered_count / total_points) * 100:.2f}%"
    }
    return stats

def create_publication_plot(plot_data: pd.DataFrame, output_file: str):
    """Creates and saves a publication-quality scatter plot."""
    if plot_data.empty:
        print("No data to plot.", file=sys.stderr)
        return

    # 1. 테마 및 컨텍스트 설정 (논문 스타일에 맞게)
    sns.set_theme(style="ticks") # palette는 scatterplot에서 직접 지정
    sns.set_context("paper", font_scale=1.5)

    fig, ax = plt.subplots(figsize=(8, 8))

    # 2. 산점도(Scatter Plot) 생성
    sns.scatterplot(
        data=plot_data,
        x="CTRL",
        y="MASK",
        hue="filtered",
        hue_order=[True, False],
        # [수정된 부분] set 대신 dictionary를 사용하여 값을 색상에 명시적으로 매핑
        palette={True: "tab:blue", False: "#BBBBBB"},
        s=50,
        alpha=0.8,
        edgecolor="w",
        linewidth=0.5,
        ax=ax
    )

    # 3. 보조선 및 강조 영역 수정
    line_color = '#888888'
    ax.axhline(THRESHOLD, color=line_color, linestyle='--', linewidth=1)
    ax.axvline(THRESHOLD, color=line_color, linestyle='--', linewidth=1)
    ax.fill_between([0, THRESHOLD], 0, THRESHOLD, color='whitesmoke', alpha=0.7, zorder=0)

    # 4. 레이블과 제목을 전문적으로 수정
    fig.suptitle("MPSE Comparison: Before and After Masking", fontsize='x-large', fontweight='bold')
    ax.set_title(f"Points retained (MPSE < {THRESHOLD}) are shown in gray", fontsize='medium', pad=10)
    ax.set_xlabel("MPSE (Before Masking)")
    ax.set_ylabel("MPSE (After Masking)")
    
    handles, labels = ax.get_legend_handles_labels()
    new_labels = [f'Filtered (MPSE > {THRESHOLD})', f'Retained (MPSE < {THRESHOLD})']
    ax.legend(handles, new_labels, title="Status", frameon=False)

    # 5. 축 범위와 테두리 정리
    ax.set_xlim(0, 1.01)
    ax.set_ylim(0, 1.01)
    ax.set_aspect('equal', adjustable='box')
    sns.despine(trim=True)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    # 6. 고해상도 벡터 파일로 저장
    output_path = Path(output_file)
    print(f"Saving plots to {output_path.parent}/")
    
    plt.savefig(output_path.with_suffix('.pdf'), bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
    plt.savefig(output_path.with_suffix('.png'), dpi=300, bbox_inches='tight')
    
    plt.close(fig)


@click.command()
@click.option('--input-file', required=True, type=click.Path(exists=True, path_type=Path), help='Input TSV file generated from analyze_bpp.')
@click.option('--output-file', required=True, type=click.Path(path_type=Path), help='Path for the output plot file.')
def main(input_file: Path, output_file: Path) -> None:
    """Main function to generate a comprehensive plot showing all data states."""
    
    print(f"Reading records from {input_file.name}...")
    try:
        data = pd.read_csv(input_file, sep='\t')
        print(f"Found {len(data)} records.")
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        sys.exit(1)

    pivoted_data = load_and_prepare(data)

    if pivoted_data.empty:
        print("Error: No data to plot after pivoting. Check input file format.")
        sys.exit(1)

    # --- 1. [수정] 데이터 포인트의 상태를 정의하는 'status' 컬럼 생성 ---
    # 기존의 success_condition 대신 모든 데이터의 상태를 구분합니다.
    conditions = [
        (pivoted_data['CTRL'] <= 0.7) & (pivoted_data['MASK'] >= 0.7), # 성공 (Success)
        (pivoted_data['CTRL'] > 0.7) | (pivoted_data['MASK'] < 0.7)  # 실패 (Failure)
    ]
    choices = ['Success', 'Failure']
    pivoted_data['status'] = pd.np.select(conditions, choices, default='Failure')

    # --- 2. [수정] 통계 계산 로직 변경 ---
    total_points = len(pivoted_data)
    success_points = len(pivoted_data[pivoted_data['status'] == 'Success'])
    failure_points = total_points - success_points
    success_rate = (success_points / total_points) * 100 if total_points > 0 else 0

    print("\n--- Data Statistics ---")
    print(f"Total points: {total_points}")
    print(f"Success points (CTRL<=0.7 & MASK>=0.7): {success_points}")
    print(f"Failure points: {failure_points}")
    print(f"Success Rate: {success_rate:.2f}%")
    print("-----------------------")
    
    # --- 3. [수정] 플롯 생성 로직 ---
    # Publication-quality 스타일 적용
    sns.set_theme(style="ticks", rc={"font.family": "sans-serif"})
    sns.set_context("paper", font_scale=1.5)
    fig, ax = plt.subplots(figsize=(10, 8))

    # hue에 'status'를 사용하여 상태별로 색상 구분
    sns.scatterplot(
        data=pivoted_data,
        x="CTRL",
        y="MASK",
        hue="status",
        hue_order=['Success', 'Failure'],
        palette={'Success': '#0072B2', 'Failure': '#BBBBBB'}, # 파란색과 회색으로 지정
        s=50,
        alpha=0.8,
        edgecolor="w",
        linewidth=0.5,
        ax=ax
    )
    
    # --- 4. [수정] 시각적 요소 정리 ---
    line_color = '#888888'
    ax.axhline(0.7, color=line_color, linestyle='--', linewidth=1)
    ax.axvline(0.7, color=line_color, linestyle='--', linewidth=1)
    
    # 성공 영역 강조
    ax.fill_between([0, 0.7], 0.7, 1, color='#0072B2', alpha=0.1, zorder=0)

    # 제목 및 레이블을 더 명확하게 수정
    ax.set_title("MPSE Comparison: Before vs. After Masking", fontsize=18, pad=15, weight='bold')
    ax.set_xlabel("MPSE (Before Masking)", fontsize=14)
    ax.set_ylabel("MPSE (After Masking)", fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)

    # 축 범위와 비율 설정
    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_aspect('equal', adjustable='box')
    
    # 범례 수정
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, title='Status', frameon=False, loc='upper left')

    # 불필요한 테두리 제거
    sns.despine(trim=True)
    plt.tight_layout()

    # --- 5. 저장 ---
    output_file.parent.mkdir(parents=True, exist_ok=True)
    # 여러 포맷으로 저장
    plt.savefig(output_file.with_suffix('.png'), dpi=300)
    plt.savefig(output_file.with_suffix('.pdf'))
    print(f"\nPlot saved to: {output_file.with_suffix('.png')} (and .pdf)")


if __name__ == '__main__':
    main()