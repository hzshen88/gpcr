# GPCR pivot depth: packing-density and WCN reproduction package

This directory reproduces the pivot-depth computation described in the
Methods section of "GPCR mechanical pivots are packed closer to the
membrane center than pivots in other polytopic membrane protein folds."

## Contents

- `data/structure_manifest.csv` — the 88-structure primary-analysis list
  (PDB ID, fold family, receptor name, UniProt accession, GPCR/non-GPCR
  group) used in the manuscript. 9 additional structures whose membrane
  orientation could not be confirmed are documented but excluded from the
  primary analysis (see manuscript Methods, "Structure set").
- `data/validation_manifest_13.csv` — a 13-structure spot-check subset
  covering all 9 fold families (plus the two structures used in the
  manuscript's Figure 1 methods-schematic panels) used by `VALIDATION.md`.
- `config/params.yaml` — depth-binning, packing-metric, orientation-tier,
  and statistical parameters, matching the manuscript Methods exactly.
- `scripts/01_fetch_structures.py` — downloads each structure's CIF file
  from RCSB PDB.
- `scripts/02_orient_and_compute_pivot.py` — orients each structure along
  the membrane normal (using the manuscript's actual 4-tier procedure --
  see below), then computes pivot depth under both local packing
  definitions (coordination-number packing, all-atom WCN).
- `results/` — script outputs land here.
- `VALIDATION.md` — a direct, structure-by-structure comparison of this
  package's output against the manuscript's published per-structure
  values, including the two explained exceptions found during validation.

## Orientation procedure

Four tiers, applied in order of preference (see manuscript Methods,
"Pivot depth estimation" for full detail):

1. **Topology-verified PCA** (primary method, 59/97 structures). Estimate
   the membrane normal via PCA (largest-variance axis) of TM-region Cα
   coordinates (TM ranges from the UniProt REST API), then verify/correct
   its sign using UniProt "Topological domain" annotations: require the
   mean depth of residues annotated Extracellular to exceed the mean depth
   of residues annotated Cytoplasmic (≥3 residues mapped to each class).
2. **PDBTM fallback** (28/97). For structures lacking usable topology
   annotations, use the precomputed membrane-boundary transformation
   matrix from the PDBTM database (TMDET algorithm) via its public XML
   API.
3. **Unverified PCA fallback** (7/97). PCA normal retained without an
   independent sign check.
4. **N/C-terminus fallback** (3/97). Coarsest tier: sign check from known
   N/C-terminus cytoplasmic localization, supplied by the caller (this
   script does not guess it silently).

**We do NOT use the OPM/PPM database anywhere in this pipeline.** An
earlier internal draft of the manuscript incorrectly described the
orientation method as OPM/PPM-based; that description has been corrected
in the manuscript (see main text, Methods, and Data availability).

## Pivot-depth computation

Both packing definitions follow the same two-step structure: (1) compute
a per-residue local-packing value using the FULL chain (all Cα atoms, or
all heavy atoms) as the neighbor pool, so no genuine spatial neighbor is
missed for residues near the TM-region boundary; (2) restrict to
TM-flagged residues only and bin by depth into a FIXED `[-21, +21]` Å
range (14 bins of width 3 Å, ≥3 residues per bin) to find the pivot (the
bin of maximum mean packing value).

- **Coordination-number packing**: per-Cα-atom count of neighboring Cα
  atoms within 10 Å (strict `<` cutoff).
- **All-atom WCN** (Weighted Contact Number): per-heavy-atom
  Σ 1/d² over all other heavy atoms within a 15 Å capture radius,
  aggregated to a per-residue mean.

Crystallization-fusion partners (T4-lysozyme, BRIL, nanobodies, etc.) are
**not** excluded by any resnum-based heuristic -- an earlier draft tried
that and it silently truncated genuine native residues for receptors
whose own numbering legitimately exceeds 1000 (e.g. 7SF7/ADGRL3). Fusion
partners are excluded implicitly and correctly instead: they fall outside
the UniProt-derived TM ranges used to build the depth-binned profile, so
they never enter it regardless of how they happen to be numbered.

This matches the manuscript's primary analysis -- see `VALIDATION.md` for
a 13-structure spot-check (covering all 9 fold families) confirming 12/13
exact matches on coordination-number packing and 11/13 on WCN, with the
one remaining discrepancy traced to a specific, understood cause (not an
unexplained bug).

## Usage

```bash
pip install gemmi numpy pandas scipy pyyaml biopython
python scripts/01_fetch_structures.py
python scripts/02_orient_and_compute_pivot.py --manifest data/structure_manifest.csv --out results/pivot_depths.csv
```

To reproduce the 13-structure validation check in `VALIDATION.md`:

```bash
python scripts/02_orient_and_compute_pivot.py --manifest data/validation_manifest_13.csv --out results/validation_test_pivot_13.csv
```

## Data availability

Per-structure pivot depth and clarity values as reported in the
manuscript are provided in the manuscript's Supplementary Table S1
(n=88, primary analysis) and can be independently reproduced with this
package (see `VALIDATION.md`).
