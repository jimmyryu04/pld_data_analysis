#!/usr/bin/env python3
"""
Scatter plot: MPSE-utp (CTRL vs MASK) for random sampling data
──────────────────────────────────────────────────────────────
python scatter_plot_mpse_random_sampling.py --input-files <input.tsv> --output-file <out.html>

• X: MPSE-utp (CTRL)
• Y: MPSE-utp (MASK)
• 색: 3utr-fwd
• hover: group 이름 (샘플 번호 포함)
"""
import sys
from pathlib import Path
import pandas as pd  # <--- 이 줄을 추가하여 문제를 해결합니다.
import plotly.express as px
import re
import click

METRIC = "MPSE-utp"   # 분석할 메트릭 이름

def load_and_prepare(df: pd.DataFrame) -> pd.DataFrame:
    """Prepares data for plotting, handling random sampling names."""
    
    # 타입 추출: 이름에서 '_CTRL' 또는 '_MASK' 부분을 찾음
    df['type'] = df['name'].str.extract(r'_(CTRL|MASK)')

    # 그룹 생성: 이름에서 '_CTRL' 또는 '_MASK' 부분을 제거하여 고유한 그룹 ID를 만듬
    df['group'] = df['name'].str.replace(r'_(CTRL|MASK)', '', regex=True)

    # 필요한 타입(CTRL, MASK)만 필터링하고, 누락된 값이 없는 행만 선택
    sub = df[df["type"].isin(["CTRL", "MASK"])].dropna(subset=['type', 'group', METRIC])
    
    if sub.empty:
        return pd.DataFrame()

    # 그룹 × 타입 → METRIC 값으로 피벗 테이블 생성
    # 중복된 인덱스가 없는지 확인 후 피벗
    sub = sub.drop_duplicates(subset=['group', 'type'])
    pvt = (sub
           .pivot(index="group",
                  columns="type",
                  values=METRIC)
           .dropna())

    # 색상 지정을 위한 3'UTR 정보 추출 (CTRL 행 기준)
    utr_map_df = df[df["type"] == "CTRL"].drop_duplicates(subset=['group'])
    utr_map = utr_map_df.set_index("group")["3utr-fwd"]
    
    pvt = pvt.merge(utr_map.rename("3utr-fwd"), left_index=True, right_index=True, how="left")
    pvt["hover"] = pvt.index  # hover 텍스트용 그룹 이름

    return pvt.reset_index()

@click.command()
@click.option('--input-files', required=True, type=click.Path(exists=True, path_type=Path), help='Input TSV file generated from analyze_bpp.')
@click.option('--output-file', required=True, type=click.Path(path_type=Path), help='Path for the output HTML plot.')
def main(input_files: Path, output_file: Path) -> None:
    """Main function to generate and save the plot."""
    
    print(f"Reading records from {input_files.name}...")
    try:
        data = pd.read_csv(input_files, sep='\t')
        print(f"Found {len(data)} records.")
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        sys.exit(1)

    plot_data = load_and_prepare(data)

    if plot_data.empty:
        print("Error: No data to plot after processing. Please check the input file content and format.")
        sys.exit(1)

    fig = px.scatter(
        plot_data,
        x="CTRL",
        y="MASK",
        color="3utr-fwd",
        hover_name="hover",
        labels={
            "CTRL":  f"MPSE (before masking)",
            "MASK":  f"MPSE (after masking)",
            "3utr-fwd": "3′UTR"
        },
        title=f"NLucP MPSE per Sample (before masking vs after masking)",
        height=700,
    )

    # 대각선 y = x 추가
    min_val = min(plot_data["CTRL"].min(), plot_data["MASK"].min()) * 0.95
    max_val = max(plot_data["CTRL"].max(), plot_data["MASK"].max()) * 1.05
    fig.add_shape(
        type="line",
        x0=min_val, y0=min_val,
        x1=max_val, y1=max_val,
        line=dict(width=1, dash="dash", color="grey"),
    )
    
    fig.update_layout(
        xaxis_title=f"MPSE (before masking)",
        yaxis_title=f"MPSE (after masking)",
        legend_title_text='3\'UTR fwd'
    )

    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(output_file)
    print(f"Plot saved to: {output_file}")


if __name__ == '__main__':
    main()