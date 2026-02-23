#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


MAJOR_NC = {"GC-AG", "AT-AC"}
CANONICAL = "GT-AG"


def load_counts(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)
    df.columns = df.columns.astype(str).str.strip()
    if "Dimer" not in df.columns:
        raise ValueError("Counts table must contain a 'Dimer' column.")
    return df.set_index("Dimer").fillna(0)


def load_totals(path: Path | None, counts: pd.DataFrame) -> pd.Series:
    if path is None:
        return counts.sum(axis=0)

    df = pd.read_csv(path)
    df.columns = df.columns.astype(str).str.strip()
    if not {"species", "total_introns"}.issubset(df.columns):
        raise ValueError("Totals CSV must have columns: species,total_introns")
    totals = pd.Series(df["total_introns"].values, index=df["species"].astype(str))
    return totals


def summarize(counts: pd.DataFrame, totals: pd.Series, mode: str) -> pd.DataFrame:
    # mode: "all" or "minor"
    counts = counts.copy()
    counts.columns = counts.columns.astype(str)

    totals = totals.reindex(counts.columns)
    if totals.isna().any():
        missing = totals[totals.isna()].index.tolist()
        raise ValueError(f"Totals missing for species: {missing}")

    if mode == "all":
        exclude = {CANONICAL}
        label = "noncanonical_all"
    elif mode == "minor":
        exclude = {CANONICAL} | MAJOR_NC
        label = "noncanonical_minor"
    else:
        raise ValueError("mode must be 'all' or 'minor'")

    dimers = [d for d in counts.index.astype(str) if d not in exclude]
    sub = counts.loc[dimers] if dimers else counts.iloc[0:0]

    occ = sub.sum(axis=0)
    types = (sub > 0).sum(axis=0)
    frac = (occ / totals).replace([np.inf, -np.inf], np.nan).fillna(0)

    canonical_counts = counts.loc[CANONICAL] if CANONICAL in counts.index else pd.Series(0, index=counts.columns)

    out = pd.DataFrame({
        "species": counts.columns,
        "total_introns": totals.values,
        "canonical_GT-AG": canonical_counts.reindex(counts.columns).values,
        f"{label}_occurrences": occ.reindex(counts.columns).values,
        f"{label}_types": types.reindex(counts.columns).values,
        f"{label}_fraction": frac.reindex(counts.columns).values,
    })
    return out


def main():
    p = argparse.ArgumentParser(
        description="Compute per-species non-canonical splice-site summary (S9-like)."
    )
    p.add_argument("--counts", required=True, help="Merged counts table (.xlsx or .csv)")
    p.add_argument("--totals", default=None, help="Totals CSV with columns species,total_introns")
    p.add_argument("--out-dir", default="outputs", help="Output folder")
    args = p.parse_args()

    counts_path = Path(args.counts)
    totals_path = Path(args.totals) if args.totals else None
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    counts = load_counts(counts_path)
    totals = load_totals(totals_path, counts)

    all_df = summarize(counts, totals, mode="all")
    minor_df = summarize(counts, totals, mode="minor")

    out_all = out_dir / "S9_noncanonical_all_per_species_summary.csv"
    out_minor = out_dir / "S9_noncanonical_minor_per_species_summary.csv"
    all_df.to_csv(out_all, index=False)
    minor_df.to_csv(out_minor, index=False)

    if totals_path is None:
        print("[WARN] totals not provided; used column sums as totals.")

    print(f"[OK] wrote: {out_all}")
    print(f"[OK] wrote: {out_minor}")


if __name__ == "__main__":
    main()
