# Provenance

## Origin of the material

1. **Operator definition.** Restated in THEOREM.md §1 from the owner's internal design document for the discretization (frozen 2026-09-19; its SHA-256 is `d894de10b356333b469bfc63c653c60bd359ec…`, recorded here only so the owner can match it to internal records). The restatement is complete for the two charts involved; no internal file is included.
2. **Derivation and proof-session scripts** (`verification/proof_session/d1, d2, d3, d5`). Produced 2026-10-31 (16:11–16:29 UTC) by an AI coding agent (Claude Opus 5.5, Claude Code, single fresh session) working from the design text only, with no import or execution of any implementation code. The scripts are included verbatim; the session's write-up was adapted into THEOREM.md and PROOF.md §1–5, §7. The session's sealed package has manifest SHA-256 `ffa33b9b9ad13e176ba1fa0ba2a3f366d9b4f573f23a7031e79186aca0c43a28` (155 files) and archive SHA-256 `4cf7414d8d9f1eef417e94643ea0cd644c96e8d24f57a1681e8215ed00e48955`.
3. **Independent check and second derivation** (`verification/independent/`, PROOF.md §6, the corollary, PRIOR_ART.md, LIMITATIONS.md). Produced 2026-10-31 (18:27–≈20:50 UTC) by a different model (Claude Fable 5.1) in a separate session that received the sealed package above, re-derived the seam row by symbolic differentiation of the assembled energy, derived the interior density and the boundary reduction, and ran a fresh literature search. The review bundle has archive SHA-256 `4d86ce7b00b15d9e1ffbef2407727390494dd29030ea0b4eafd60da473e05114`.
4. **Owner.** Michael Blumhardt directed both sessions, holds the internal records, and is responsible for this release.

## What is deliberately omitted

- Implementation code of the discretization and of the instrument it was a candidate for.
- The owner's numerical result records (regression coefficients, band-functional measurements). They were consistent with the exact results (LIMITATIONS.md) but are not proof steps and are not released.
- Internal custody paths, session identifiers and governance records.

## Reproducibility

`verification/run_all.sh` regenerates every symbolic result from this repository alone; expected outputs are in `verification/expected_output/`. The independent scripts were run from a clean copy of this tree before release and reproduced the proof-session outputs byte-for-byte apart from timings.

## Revision history

- r1 (2026-10-31): initial draft.
- r2: adjoint-defect correction; conditional general form, examples and counterexamples; Track A observation; Linders 2018 examined.
- r3: residual-norm corollary split into conservative (single-witness) and symmetric (pair) bounds with |W′|-weighted norms; the symmetric part of the general form restated as a pair bound; single-witness symmetric cancellation remark added; full verification suite rerun (1 m 33 s) with outputs matching `verification/expected_output/`.
- r4: editorial only — PRIOR_ART.md opening paragraph now says one unexamined text, consistent with §4 and the README; full-suite run log of the r3 rerun added as `verification/RUN_LOG_2026-10-31.txt`.

## Manifest

`MANIFEST.sha256` lists the SHA-256 of every file in the repository except itself. Regenerate with:

```
find . -type f ! -path './.git/*' ! -name MANIFEST.sha256 -print0 | sort -z | xargs -0 sha256sum > provenance/MANIFEST.sha256
```

## Public-release edits (v1.0.0)

- `verification/proof_session/d1_seam_series.py`: SHA-256 before `6a399a7b7d9078333ac3d633218bdbd454814dd78a1036d2a2b9b70284eaae5c`, after `a7239de25cbbb04c9a9205039029b6fa9979ecf326df7bec99308fdf1f915171`. Docstring wording edited for public release; no code changed.
- `verification/proof_session/d2_fermi_reduction.py`: SHA-256 before `a1da003e733927ed8daa175556ee2e5d5d28a363dd89ef718b28768eaedf04a9`, after `12698d018fe856e4dbf4e2fb4397862b1c95a8f2b138e6ca5b803d3846efa8d6`. Docstring wording edited for public release; no code changed.
- `verification/independent/independent_seam.py`: SHA-256 before `e056ceeef76f8ca64dcbe45167e4df844b83c70c1b0321bd656e4d3ae95950dc`, after `b6ef342aabcbc200f9a7adc2dc53b3a9c15400dbfe73f59580f081c4c8e0e0ba`. Docstring wording edited for public release; no code changed.
- `verification/proof_session/d5_witness_bounds.py`: SHA-256 before `c51847efabed9586828bba9b1e2e7f2da8ac22638114e1561b2a69abf5f1278c`, after `2e79d0675305685dac6962b4d4c07f2b9893a1671bbfb66d317d03cfa7072ba6`. Docstring wording edited for public release; no code changed.
- Added `LICENSE` (CC BY 4.0, for all Markdown files), `verification/LICENSE` (MIT, for the code in `verification/`) and `CITATION.cff`.
- `README.md`: the License section now states the licences and points to the two LICENSE files; the Citation section now gives the repository URL and tag v1.0.0 in place of the publication placeholder.
- Verification rerun after these edits: `verification/run_all.sh` in a clean virtual environment (Python 3.13.11, sympy 1.14.0, mpmath 1.3.0, numpy 2.5.4; x86_64 Linux), wall time 38.1 s; all six outputs matched expected_output/ apart from timing lines.
- Manifest regenerated after these edits; no manifest-covered file was changed afterwards.
