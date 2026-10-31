#!/usr/bin/env bash
# Reproduce every symbolic check from a clean environment.
# Requirements: python3 (>=3.10), sympy, mpmath.  numpy is used only by independent_seam.py for nothing heavy.
# All scripts are exact-arithmetic (sympy rationals/surds); total run time ~2-3 min on a laptop.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output
echo "== proof-session checks (verbatim from the proof package) =="
python3 -I proof_session/d1_seam_series.py        output/d1_seam_series.json      | tee output/d1_seam_series.txt
python3 -I proof_session/d2_fermi_reduction.py    output/d2_fermi_reduction.json  proof_session | tee output/d2_fermi_reduction.txt
python3 -I proof_session/d3_interior_series.py    | tee output/d3_interior_series.txt
python3 -I proof_session/d5_witness_bounds.py     | tee output/d5_witness_bounds.txt
echo "== independent checks (review session, different derivation routes) =="
python3 -I independent/independent_seam.py        | tee output/independent_seam.txt
python3 -I independent/independent_boundary.py    | tee output/independent_boundary.txt
echo
echo "Compare output/*.txt with expected_output/*.txt (timings excepted)."
