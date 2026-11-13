# Old, unused code

Kept for the record. Nothing here is used by the paper's benchmarks or figures,
and paths inside these files still point at their old locations, so they will
not run from here as is.

- `workflows/create/`: the first creation benchmark (one Snakefile over sample
  counts and thread counts), replaced by `workflows/cohort/` and
  `workflows/single/`. Its `timed_trials.py` and `bigwig_scan.py` moved to
  `workflows/lib/`.
- over10 (gnomAD-style fraction of samples at >= 10x): `workflows/create/over10_*.py`,
  `workflows/create/proposed_api/`, `workflows/cohort/over10_d4.py`,
  `workflows/cohort/summary_*.py` and `config/cohort-bigwig-over10.yaml`.
- `config/1kg-*.yaml`: configs for the first creation benchmark.
- `tests/`: tests for the code above.
