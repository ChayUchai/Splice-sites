#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


BASES = ["A", "C", "G", "T"]


def is_acgt_combo(combo: str) -> bool:
    s = str(combo).replace("-", "")
    return len(s) == 4 and set(s) <= set("ACGT")


def combo_base_freqs(combo: str) -> dict[str, float]:
    s = str(combo).replace("-", "")
    freq = {b: 0.0 for b in BASES}
    for ch in s:
        if ch in freq:
            freq[ch] += 1.0
    for b in BASES:
        freq[b] /= 4.0
    return freq


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


def long_short_counts(df: pd.DataFrame, long_thr: int, pseudo: float) -> pd.DataFrame:
    df = df.copy()
    df = df[df["combo"].apply(is_acgt_combo)]
    df["intron_len"] = pd.to_numeric(df["intron_len"], errors="coerce")
    df = df.dropna(subset=["intron_len"])

    df["is_long"] = df["intron_len"] >= long_thr

    grp = df.groupby(["organism", "combo", "is_long"]).size().unstack(fill_value=0)
    grp.columns = ["short", "long"] if False in grp.columns and True in grp.columns else [str(c) for c in grp.columns]

    if "short" not in grp.columns:
        grp["short"] = 0
    if "long" not in grp.columns:
        grp["long"] = 0

    grp = grp.reset_index()
    grp["ratio_long_short"] = (grp["long"] + pseudo) / (grp["short"] + pseudo)
    return grp


def rank_and_spearman(one: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # one: combo-level table for a single organism (or ALL_POOLED)
    t = one.copy()
    t = t.sort_values("ratio_long_short", ascending=False).reset_index(drop=True)
    t["rank"] = np.arange(1, len(t) + 1)

    # base freqs
    bf = t["combo"].apply(combo_base_freqs).apply(pd.Series)
    bf.columns = [f"freq_{c}" for c in bf.columns]
    t = pd.concat([t, bf], axis=1)

    rows = []
    for b in BASES:
        r, p = spearmanr(t["rank"], t[f"freq_{b}"])
        rows.append({"base": b, "spearman_r": r, "p_value": p})
    stats = pd.DataFrame(rows)
    return t, stats


def main():
    p = argparse.ArgumentParser(
        description="Long vs short introns per splice-site combo + base-frequency correlations."
    )
    p.add_argument("--input-dir", required=True, help="Folder with intron length CSV files")
    p.add_argument("--pattern", default="*.csv", help="Glob pattern (default: *.csv)")
    p.add_argument("--long-thr", type=int, default=5000, help="Threshold for long introns (bp)")
    p.add_argument("--pseudo", type=float, default=1.0, help="Pseudocount for ratio")
    p.add_argument("--out-dir", default="outputs", help="Output folder")
    args = p.parse_args()

    input_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_intron_tables(input_dir, args.pattern)

    # per organism
    ls = long_short_counts(df, long_thr=args.long_thr, pseudo=args.pseudo)

    ranked_all = []
    stats_all = []

    for org, sub in ls.groupby("organism"):
        ranked, stats = rank_and_spearman(sub[["combo", "short", "long", "ratio_long_short"]])
        ranked.insert(0, "organism", org)
        stats.insert(0, "organism", org)
        ranked_all.append(ranked)
        stats_all.append(stats)

    # pooled (ALL_POOLED)
    pooled = df.copy()
    pooled["organism"] = "ALL_POOLED"
    pooled_ls = long_short_counts(pooled, long_thr=args.long_thr, pseudo=args.pseudo)
    pooled_ranked, pooled_stats = rank_and_spearman(
        pooled_ls[["combo", "short", "long", "ratio_long_short"]]
    )
    pooled_ranked.insert(0, "organism", "ALL_POOLED")
    pooled_stats.insert(0, "organism", "ALL_POOLED")
    ranked_all.insert(0, pooled_ranked)
    stats_all.insert(0, pooled_stats)

    ranked_df = pd.concat(ranked_all, ignore_index=True)
    stats_df = pd.concat(stats_all, ignore_index=True)

    out_ranked = out_dir / "all_combos_ranked_longshort.csv"
    out_stats = out_dir / "basefreq_spearman.csv"
    ranked_df.to_csv(out_ranked, index=False)
    stats_df.to_csv(out_stats, index=False)

    print(f"[OK] wrote: {out_ranked}")
    print(f"[OK] wrote: {out_stats}")


if __name__ == "__main__":
    main()
