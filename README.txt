Splice-site analysis scripts

Python scripts for analysis of splice-site combinations in plant
genomes. The repository includes tools for merging per-species counts,
correlation analysis, non-canonical splice-site summary, intron length
distribution, and ranking of combinations.

  ---------------------
  Structure
  ---------------------
  Requirements

  Python 3.9+
  recommended.

  Install packages:

  pip install -r
  requirements.txt
  ---------------------

Input data

1.  Per-species splice-site counts

CSV format:

Dimer,count GT-AG,12345 GC-AG,67 AT-AC,5

Each file corresponds to one species.

2.  Total intron numbers

Used for normalization. Species names must match column names in merged
tables.

3.  Intron length tables

CSV files per species:

length,combo,combo_total

Combo format: NN-NN (example: GT-AG).

  ---------------------
  Workflow
  ---------------------
  Notes

  Use relative paths or
  CLI arguments for
  input/output. Do not
  store large raw data
  in the repository.
  Keep only small
  example datasets in
  examples/.
  ---------------------

License

TBD
