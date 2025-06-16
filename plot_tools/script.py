#!/usr/bin/env python3
"""
CTRL vs MASK metric counter (+ totals)
──────────────────────────────────────
python count_ctrl_vs_mask.py <input.tsv> <output.tsv>
"""

import sys
import pandas as pd

def main(in_tsv: str, out_tsv: str) -> None:
    df = pd.read_csv(in_tsv, sep="\t")

    need = ["group", "type", "3utr-fwd", "bp-utp-ms", "MPSE-utp"]
    df = df[need]

    rows = []
    for utr, sub in df.groupby("3utr-fwd"):
        pvt = sub.pivot(index="group",
                        columns="type",
                        values=["bp-utp-ms", "MPSE-utp"])

        bp_tbl = pvt["bp-utp-ms"][["CTRL", "MASK"]].dropna()
        mp_tbl = pvt["MPSE-utp"][["CTRL", "MASK"]].dropna()

        # 두 metric 모두에 대해 CTRL·MASK가 존재하는 그룹만 세고 싶다면:
        # valid_idx = bp_tbl.index.intersection(mp_tbl.index)
        # bp_tbl, mp_tbl = bp_tbl.loc[valid_idx], mp_tbl.loc[valid_idx]

        rows.append({
            "3utr-fwd": utr,
            "total_groups":              bp_tbl.shape[0],          # 비교에 사용된 그룹 수
            "bp_utp_ms_CTRL_lt_MASK":   (bp_tbl["CTRL"] <  bp_tbl["MASK"]).sum(),
            "bp_utp_ms_equal":          (bp_tbl["CTRL"] == bp_tbl["MASK"]).sum(),
            "MPSE_utp_CTRL_lt_MASK":    (mp_tbl["CTRL"] <  mp_tbl["MASK"]).sum(),
            "MPSE_utp_equal":           (mp_tbl["CTRL"] == mp_tbl["MASK"]).sum(),
        })

    out_df = pd.DataFrame(rows)

    # ──────────────── 전체 합계(TOTAL 행) 추가 ──────────────── #
    total_row = {
        "3utr-fwd":      "TOTAL",
        "total_groups":  out_df["total_groups"].sum(),
        "bp_utp_ms_CTRL_lt_MASK": out_df["bp_utp_ms_CTRL_lt_MASK"].sum(),
        "bp_utp_ms_equal":        out_df["bp_utp_ms_equal"].sum(),
        "MPSE_utp_CTRL_lt_MASK":  out_df["MPSE_utp_CTRL_lt_MASK"].sum(),
        "MPSE_utp_equal":         out_df["MPSE_utp_equal"].sum(),
    }
    out_df = pd.concat([out_df, pd.DataFrame([total_row])], ignore_index=True)

    out_df.to_csv(out_tsv, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python count_ctrl_vs_mask.py <input.tsv> <output.tsv>")
    main(sys.argv[1], sys.argv[2])