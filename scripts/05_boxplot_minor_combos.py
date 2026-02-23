#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_ranked(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = df.columns.astype(str).str.strip()
    if "organism" not in df.columns or "combo" not in df.columns:
        raise ValueError("Ranked CSV must have columns: organism, combo")
    if "ratio_long_short" not in df.columns and "ratio" in df.columns:
        df["ratio_long_short"] = df["ratio"]
    return df


def main():
    p = argparse.ArgumentParser(
        description="Boxplot of counts for selected minor splice-site combos (log scale)."
    )
    p.add_argument("--ranked", required=True, help="CSV from 04_long_short_rank_basefreq.py")
    p.add_argument("--organism", default="ALL_POOLED", help="Which organism to plot (default: ALL_POOLED)")
    p.add_argument("--focus", default="AT-TT,GT-TT,CT-AC,CT-GC,GC-TA",
                   help="Comma-separated list of combos to highlight")
    p.add_argument("--value-col", default="combo_total",
                   help="Column to plot if present; fallback to (short+long)")
    p.add_argument("--out-png", default="outputs/splice_site_boxplot_counts_from_minor_log.png",
                   help="Output PNG path")
    args = p.parse_args()

    ranked = load_ranked(Path(args.ranked))
    sub = ranked[ranked["organism"].astype(str) == args.organism].copy()
    if sub.empty:
        raise ValueError(f"No rows for organism={args.organism}")

    focus = [x.strip() for x in str(args.focus).split(",") if x.strip()]
    sub["category"] = np.where(sub["combo"].isin(focus), sub["combo"], "others")

    if args.value_col in sub.columns:
        sub["value"] = pd.to_numeric(sub[args.value_col], errors="coerce")
    else:
        # fallback: build a count proxy
        if {"short", "long"}.issubset(sub.columns):
            sub["value"] = pd.to_numeric(sub["short"], errors="coerce").fillna(0) + pd.to_numeric(sub["long"], errors="coerce").fillna(0)
        else:
            raise ValueError(f"Column {args.value_col} not found and no short/long columns available")

    sub = sub.dropna(subset=["value"])
    sub = sub[sub["value"] > 0]

    order = focus + (["others"] if "others" in sub["category"].unique() else [])
    data = [sub.loc[sub["category"] == cat, "value"].values for cat in order]

    plt.figure(figsize=(10, 4), dpi=200)
    plt.boxplot(data, labels=order, showfliers=False)
    plt.yscale("log")
    plt.xlabel("splice-site combo")
    plt.ylabel("counts (log scale)")
    plt.tight_layout()

    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"[OK] wrote: {out_png}")


if __name__ == "__main__":
    main()
