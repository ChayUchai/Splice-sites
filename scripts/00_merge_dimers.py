#!/usr/bin/env python3
import argparse
from pathlib import Path

import pandas as pd


def read_dimer_csv(path: Path) -> pd.Series:
    """Return a Series indexed by Dimer with counts."""
    df = pd.read_csv(path)
    # Allow either explicit columns or first two columns
    cols = [c.strip() for c in df.columns]
    df.columns = cols
    if "Dimer" in df.columns and "count" in df.columns:
        dimer_col, count_col = "Dimer", "count"
    else:
        dimer_col, count_col = df.columns[0], df.columns[1]
    s = df[[dimer_col, count_col]].copy()
    s[dimer_col] = s[dimer_col].astype(str).str.strip()
    s[count_col] = pd.to_numeric(s[count_col], errors="coerce").fillna(0)
    return s.set_index(dimer_col)[count_col]


def merge_folder(input_dir: Path, pattern: str) -> pd.DataFrame:
    files = sorted(input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found: {input_dir}/{pattern}")

    merged = None
    for f in files:
        species = f.stem
        series = read_dimer_csv(f).rename(species)
        if merged is None:
            merged = series.to_frame()
        else:
            merged = merged.join(series, how="outer")

    merged = merged.fillna(0)
    merged.index.name = "Dimer"
    return merged.reset_index()


def main():
    p = argparse.ArgumentParser(
        description="Merge per-species splice-site dimer count CSVs into one table."
    )
    p.add_argument("--input-dir", required=True, help="Folder with per-species CSV files")
    p.add_argument("--pattern", default="*.csv", help="Glob pattern (default: *.csv)")
    p.add_argument("--out-xlsx", default="outputs/dimers_merged.xlsx", help="Output XLSX path")
    p.add_argument("--out-csv", default="outputs/dimers_merged.csv", help="Output CSV path")
    args = p.parse_args()

    input_dir = Path(args.input_dir)
    out_xlsx = Path(args.out_xlsx)
    out_csv = Path(args.out_csv)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    merged = merge_folder(input_dir, args.pattern)
    merged.to_excel(out_xlsx, index=False)
    merged.to_csv(out_csv, index=False)

    print(f"[OK] merged files: {input_dir} ({args.pattern})")
    print(f"[OK] wrote: {out_xlsx}")
    print(f"[OK] wrote: {out_csv}")


if __name__ == "__main__":
    main()
