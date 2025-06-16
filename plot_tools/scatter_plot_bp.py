#!/usr/bin/env python3
"""
Scatter plot: bp-utp-ms (CTRL vs MASK)
──────────────────────────────────────
python plot_bp_utp_ms_scatter.py <input.tsv> [out.html]

• X: bp-utp-ms (CTRL)
• Y: bp-utp-ms (MASK)
• 색: 3utr-fwd
• hover: group 이름
"""

import sys
from pathlib import Path
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

def load_and_prepare(tsv: str) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t")

    # ───── name → group / type 컬럼 파싱(없을 경우) ───── #
    if "type" not in df.columns or "group" not in df.columns:
        # name 예시: NLucP_1E_HBA30_lam2_CTRL
        df["type"]  = df["name"].str.split("_").str[-1]          # CTRL / CFLD / MASK
        df["group"] = df["name"].str.rsplit("_", n=1).str[0]     # 나머지 = 그룹 ID

    # CTRL·MASK 두 타입에 대해서만 사용
    sub = df[df["type"].isin(["CTRL", "MASK"])]

    # 그룹 × 타입 → bp-utp-ms 피벗
    pvt = sub.pivot(index="group",
                    columns="type",
                    values="bp-utp-ms").dropna()

    # 3′UTR(색상용) 정보는 CTRL 행의 값을 따옴
    utr_map = (sub[sub["type"] == "CTRL"]
               .set_index("group")["3utr-fwd"])
    pvt["3utr-fwd"] = utr_map

    # hover 정보
    pvt["hover"] = pvt.index

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
            "CTRL": "bp-utp-ms (contol)",
            "MASK": "bp-utp-ms (mask)",
            "3utr-fwd": "3′UTR"
        },
        title="control vs mask bp-utp-ms per group",
        height=700,
    )

    # y = x 기준선 추가(시각적으로 CTRL과 MASK가 같은 지점)
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
        print(f"🔗 Plot saved 👉 {Path(out_html).resolve()}")
    else:
        fig.show()

if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit("Usage: python plot_bp_utp_ms_scatter.py <input.tsv> [output.html]")
    main(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else None)