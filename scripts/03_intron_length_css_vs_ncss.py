#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_intron_tables(input_dir: Path, pattern: str) -> pd.DataFrame:
    files = sorted(input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found: {input_dir}/{pattern}")

    dfs = []
    for f in files:
        df = pd.read_csv(f, header=None, names=["intron_len", "combo", "combo_total"])
        df["organism"] = f.stem
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True)


def is_acgt_combo(combo: str) -> bool:
    s = str(combo).replace("-", "")
    return len(s) == 4 and set(s) <= set("ACGT")


def make_curves(data: pd.DataFrame, max_len: int, bin_size: int) -> pd.DataFrame:
    data = data.copy()
    data = data[data["combo"].apply(is_acgt_combo)]
    data["class"] = np.where(data["combo"] == "GT-AG", "CSS", "NCSS")

    data["intron_len"] = pd.to_numeric(data["intron_len"], errors="coerce")
    data = data.dropna(subset=["intron_len"])
    data = data[(data["intron_len"] >= 0) & (data["intron_len"] <= max_len)]

    bins = np.arange(0, max_len + bin_size, bin_size)
    data["bin"] = pd.cut(data["intron_len"], bins=bins, include_lowest=True, right=False)

    rows = []
    for cls, sub in data.groupby("class"):
        counts = sub["bin"].value_counts().sort_index()
        total = counts.sum()
        if total == 0:
            continue
        mids = np.array([iv.left + (iv.length / 2) for iv in counts.index.categories])
        freq = (counts.values / total) * 100.0
        rows.append(pd.DataFrame({
            "class": cls,
            "bin_mid_bp": mids,
            "freq_percent": freq,
        }))

    if not rows:
        return pd.DataFrame(columns=["class", "bin_mid_bp", "freq_percent"])

    return pd.concat(rows, ignore_index=True)


def plot_curves(curves: pd.DataFrame, out_png: Path, title: str):
    plt.figure(figsize=(8, 5), dpi=200)
    for cls, sub in curves.groupby("class"):
        y = sub["freq_percent"].rolling(window=7, center=True, min_periods=1).mean()
        plt.plot(sub["bin_mid_bp"], y, linewidth=1.5, label=cls.lower())
    plt.xlabel("intron length [bp]")
    plt.ylabel("frequency [%]")
    plt.title(title)
    if len(curves["class"].unique()) > 1:
        plt.legend()
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    p = argparse.ArgumentParser(
        description="Plot intron length distribution: CSS (GT-AG) vs NCSS (others)."
    )
    p.add_argument("--input-dir", required=True, help="Folder with intron length CSV files")
    p.add_argument("--pattern", default="*.csv", help="Glob pattern (default: *.csv)")
    p.add_argument("--max-len", type=int, default=5000, help="Max intron length to include")
    p.add_argument("--bin-size", type=int, default=5, help="Bin size in bp")
    p.add_argument("--out-prefix", default="outputs/fig3_css_vs_ncss_length_distribution",
                   help="Output prefix (without extension)")
    p.add_argument("--title", default="CSS vs NCSS intron length distribution", help="Plot title")
    args = p.parse_args()

    input_dir = Path(args.input_dir)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    data = load_intron_tables(input_dir, args.pattern)
    curves = make_curves(data, max_len=args.max_len, bin_size=args.bin_size)

    out_png = out_prefix.with_suffix(".png")
    out_csv = out_prefix.with_suffix(".csv")

    plot_curves(curves, out_png, args.title)
    curves.to_csv(out_csv, index=False)

    print(f"[OK] wrote: {out_png}")
    print(f"[OK] wrote: {out_csv}")


if __name__ == "__main__":
    main()
