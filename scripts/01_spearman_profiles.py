#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt


def load_counts(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)

    df.columns = df.columns.astype(str).str.strip()
    if "Dimer" not in df.columns:
        raise ValueError("Counts table must contain a 'Dimer' column.")
    df = df.set_index("Dimer").fillna(0)
    return df


def load_totals(path: Path | None, counts: pd.DataFrame) -> pd.Series:
    if path is None:
        # Fallback: assume totals == column sums
        totals = counts.sum(axis=0)
        return totals

    df = pd.read_csv(path)
    df.columns = df.columns.astype(str).str.strip()
    if not {"species", "total_introns"}.issubset(df.columns):
        raise ValueError("Totals CSV must have columns: species,total_introns")
    totals = pd.Series(df["total_introns"].values, index=df["species"].astype(str))
    return totals


def spearman_matrix(norm: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    species = list(norm.columns)
    n = len(species)
    corr = pd.DataFrame(np.eye(n), index=species, columns=species, dtype=float)
    pval = pd.DataFrame(np.zeros((n, n)), index=species, columns=species, dtype=float)

    for i in range(n):
        for j in range(i + 1, n):
            r, p = spearmanr(norm.iloc[:, i], norm.iloc[:, j])
            corr.iat[i, j] = corr.iat[j, i] = r
            pval.iat[i, j] = pval.iat[j, i] = p

    return corr, pval


def plot_heatmap(corr: pd.DataFrame, out_png: Path, title: str):
    plt.figure(figsize=(10, 8))
    plt.imshow(corr.values, aspect="auto")
    plt.colorbar(label="Spearman r")
    plt.xticks(range(len(corr.columns)), corr.columns, rotation=90, fontsize=7)
    plt.yticks(range(len(corr.index)), corr.index, fontsize=7)
    plt.title(title)
    plt.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    p = argparse.ArgumentParser(
        description="Compute Spearman correlation between species splice-site profiles."
    )
    p.add_argument("--counts", required=True, help="Merged counts table (.xlsx or .csv)")
    p.add_argument("--totals", default=None, help="Totals CSV with columns species,total_introns")
    p.add_argument("--out-prefix", default="outputs/spearman", help="Prefix for outputs")
    p.add_argument("--title", default="Splice-site profile similarity", help="Heatmap title")
    args = p.parse_args()

    counts_path = Path(args.counts)
    totals_path = Path(args.totals) if args.totals else None
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    counts = load_counts(counts_path)
    totals = load_totals(totals_path, counts)

    # Align and normalize
    common = [c for c in counts.columns if c in totals.index]
    if not common:
        raise ValueError("No matching species names between counts columns and totals.")
    counts = counts[common]
    totals = totals[common].astype(float)

    norm = counts.div(totals, axis=1).replace([np.inf, -np.inf], np.nan).fillna(0)

    corr, pval = spearman_matrix(norm)

    corr_csv = out_prefix.parent / f"{out_prefix.name}_matrix.csv"
    pval_csv = out_prefix.parent / f"{out_prefix.name}_pvalues.csv"
    png = out_prefix.parent / f"{out_prefix.name}_heatmap.png"

    corr.to_csv(corr_csv)
    pval.to_csv(pval_csv)
    plot_heatmap(corr, png, args.title)

    if totals_path is None:
        print("[WARN] totals not provided; used column sums as totals.")

    print(f"[OK] wrote: {corr_csv}")
    print(f"[OK] wrote: {pval_csv}")
    print(f"[OK] wrote: {png}")


if __name__ == "__main__":
    main()
