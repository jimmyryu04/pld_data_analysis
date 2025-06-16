#!/usr/bin/env python3
"""
Scatter plot: MPSE-utp (CTRL vs MASK)
─────────────────────────────────────
python plot_mpse_utp_scatter.py <input.tsv> [out.html]

• X: MPSE-utp (CTRL)
• Y: MPSE-utp (MASK)
• 색: 3utr-fwd
• hover: group 이름
"""

import sys
from pathlib import Path
import pandas as pd
import plotly.express as px

METRIC = "MPSE-utp"   # 변경된 메트릭 이름만 여기서 정의

def load_and_prepare(tsv: str) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t")

    # name → group/type 컬럼이 이미 있으면 건너뜀
    if "type" not in df.columns or "group" not in df.columns:
        df["type"]  = df["name"].str.split("_").str[-1]          # CTRL / CFLD / MASK
        df["group"] = df["name"].str.rsplit("_", n=1).str[0]

    sub = df[df["type"].isin(["CTRL", "MASK"])]

    # 그룹 × 타입 → MPSE-utp 피벗 테이블
    pvt = (sub
           .pivot(index="group",
                  columns="type",
                  values=METRIC)
           .dropna())

    # 색상용 3′UTR 정보(CTRL 행 기준)
    utr_map = (sub[sub["type"] == "CTRL"]
               .set_index("group")["3utr-fwd"])
    pvt["3utr-fwd"] = utr_map
    pvt["hover"] = pvt.index               # hover용 그룹 이름

    return pvt.reset_index()

def main(in_tsv: str, out_html: str | None = None) -> None:
    data = load_and_prepare(in_tsv)

    fig = px.scatter(
        data,
        x="CTRL",
        y="MASK",
        color="3utr-fwd",
        hover_name="hover",
        labels={
            "CTRL":  f"{METRIC} (control)",
            "MASK":  f"{METRIC} (mask)",
            "3utr-fwd": "3′UTR"
        },
        title="NLucP MPSE-uridine per group",
        height=700,
    )

    # 대각선 y = x
    fig.add_shape(
        type="line",
        x0=data["CTRL"].min(), y0=data["CTRL"].min(),
        x1=data["CTRL"].max(), y1=data["CTRL"].max(),
        line=dict(width=1, dash="dash"),
        layer="below"
    )

    fig.update_layout(legend_title_text="3′UTR-fwd")

    if out_html:
        fig.write_html(out_html)
        print(f"Plot saved to {Path(out_html).resolve()}")
    else:
        fig.show()

if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit("Usage: python plot_mpse_utp_scatter.py <input.tsv> [output.html]")
    main(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else None)