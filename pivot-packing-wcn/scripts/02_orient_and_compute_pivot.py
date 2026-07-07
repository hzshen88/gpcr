#!/usr/bin/env python3
"""
Step 2: Orient each structure along the membrane normal using the SAME
three-tier procedure used in the manuscript's primary analysis, then compute
pivot depth under two independent local packing definitions (coordination-
number packing, all-atom WCN).

Orientation procedure (applied in order of preference; see manuscript
Methods, "Pivot depth estimation"):
  (1) Topology-verified PCA. Compute a candidate membrane normal via PCA
      (largest-variance axis) of TM-region C-alpha coordinates, using
      UniProt REST-API "Transmembrane" features to identify the TM region.
      Verify/correct the sign of the normal using UniProt "Topological
      domain" annotations: require the mean depth of residues annotated
      Extracellular to exceed the mean depth of residues annotated
      Cytoplasmic (>=3 residues mapped to each class). This is the primary
      method (covers 59/97 structures in the manuscript dataset).
  (2) PDBTM fallback. For structures lacking usable topology annotations,
      use the pre-computed membrane-boundary transformation matrix from the
      PDBTM database (Protein Data Bank of Transmembrane Proteins; TMDET
      algorithm) via its public XML API. Covers 28/97 structures.
  (3) N/C-terminus fallback. For structures resolved by neither method,
      fall back to a coarser sign check based on the known cytoplasmic
      localization of the N- or C-terminus (lowest-confidence tier; 3/97
      structures in the manuscript dataset). This script emits a warning
      and requires the user to supply the terminus side via the manifest
      when this tier is reached, rather than guessing silently.

We do NOT use the OPM/PPM database anywhere in this pipeline; an earlier
internal draft of the manuscript incorrectly described the orientation
method as OPM/PPM-based; that description has been corrected (see main
manuscript, Methods).

Structures whose orientation cannot be resolved by any of the three tiers
are reported with orient_confidence="unresolved" and are excluded from the
primary analysis (see manuscript, "8 structures excluded").

Usage:
    python 02_orient_and_compute_pivot.py [--manifest PATH] [--cif-dir PATH] [--config PATH] [--out PATH]
"""
import argparse
import json
import os
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np
import pandas as pd
import yaml
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree

try:
    import gemmi
except ImportError:
    sys.exit("This script requires gemmi. Install with: pip install gemmi")

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


# ---------------------------------------------------------------------------
# UniProt / PDBTM fetch helpers
# ---------------------------------------------------------------------------

def _jget(url, headers=None, retries=3, timeout=20):
    headers = headers or {"Accept": "application/json"}
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers=headers)
            raw = urllib.request.urlopen(req, timeout=timeout).read()
            if raw:
                return raw
        except Exception:
            time.sleep(0.5)
    return None


def fetch_uniprot_tm_and_topology(acc):
    """Return (tm_ranges, topo_domains) for a UniProt accession.

    tm_ranges: list of (start, end) 1-based residue ranges annotated
        "Transmembrane".
    topo_domains: list of (start, end, description) for "Topological
        domain" features (description is typically "Extracellular" or
        "Cytoplasmic").
    """
    raw = _jget(f"https://rest.uniprot.org/uniprotkb/{acc}.json")
    if raw is None:
        return [], []
    d = json.loads(raw)
    tm_ranges, topo = [], []
    for f in d.get("features", []):
        s = f["location"]["start"]["value"]
        e = f["location"]["end"]["value"]
        if s is None or e is None:
            continue
        if f["type"] == "Transmembrane":
            tm_ranges.append((int(s), int(e)))
        elif f["type"] == "Topological domain":
            topo.append((int(s), int(e), f.get("description", "")))
    return tm_ranges, topo


def fetch_pdbtm_matrix(pdb_id):
    """Return the PDBTM TMDET membrane-boundary transformation matrix for a
    PDB entry, or None if PDBTM has no entry / no membrane block for it.
    """
    raw = _jget(
        f"https://pdbtm.unitmp.org/api/v1/entry/{pdb_id.lower()}.xml",
        headers={"User-Agent": "Mozilla/5.0"},
    )
    if raw is None:
        return None
    try:
        text = raw.decode("iso-8859-1")
        root = ET.fromstring(text)
    except Exception:
        return None
    if root.get("TMP") != "yes":
        return None
    ns = {"p": "https://pdbtm.unitmp.org"}
    mem = root.find("p:MEMBRANE", ns)
    if mem is None:
        return None
    tmatrix = mem.find("p:TMATRIX", ns)
    if tmatrix is None:
        return None
    rows = {}
    for rowname in ("ROWX", "ROWY", "ROWZ"):
        e = tmatrix.find(f"p:{rowname}", ns)
        if e is None:
            return None
        rows[rowname] = (float(e.get("X")), float(e.get("Y")), float(e.get("Z")), float(e.get("T")))
    return rows


def pdbtm_depth(coords, tmatrix):
    zx, zy, zz, zt = tmatrix["ROWZ"]
    return coords[:, 0] * zx + coords[:, 1] * zy + coords[:, 2] * zz + zt


# ---------------------------------------------------------------------------
# Structure parsing
# ---------------------------------------------------------------------------

def list_chains(cif_path):
    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    return [(c.name, len([r for r in c if r.find_atom("CA", "\0") is not None])) for c in st[0]]


def get_ca_and_heavy(cif_path, chain_id=None):
    """Extract CA and heavy-atom coordinates for one chain.

    IMPORTANT: for multi-chain depositions (e.g. a GPCR solved together with
    a bound G-protein heterotrimer and/or nanobody), `chain_id` MUST be set
    to the receptor's own chain -- pooling all chains together will include
    non-membrane-protein atoms in the packing calculation. Verify chain
    identity against the CIF's own _entity.pdbx_description /
    _entity_poly.pdbx_strand_id records, not by guessing from the chain
    letter (see VALIDATION.md for a worked example where this went wrong).
    """
    st = gemmi.read_structure(cif_path)
    st.setup_entities()
    model = st[0]
    chains = [c for c in model if chain_id is None or c.name == chain_id]
    if not chains:
        return None
    # if chain_id not specified, use the largest protein chain
    if chain_id is None:
        chains = sorted(chains, key=lambda c: len([r for r in c if r.find_atom("CA", "\0") is not None]), reverse=True)
    ch = chains[0]
    ca_coords, ca_resnums = [], []
    heavy_coords, heavy_resnums = [], []
    for res in ch:
        # Deduplicate alternate-location (altloc) conformers: keep exactly one
        # atom per atom-name, preferring blank altloc, else the first altloc
        # encountered (typically 'A') -- matches Biopython's default disordered-
        # atom selection, which the manuscript's primary analysis pipeline used.
        # gemmi's find_atom(name, "\0") silently returns None (skipping the
        # residue entirely) when EVERY conformer of that atom has a non-blank
        # altloc, which under-counts residues for structures with alternate
        # conformers (e.g. 5KUK) -- this explicit dedup avoids that.
        seen_names = {}
        for atom in res:
            if atom.name not in seen_names or atom.altloc in ("\0", "A", ""):
                if atom.name not in seen_names or seen_names[atom.name].altloc not in ("\0", "A", ""):
                    seen_names[atom.name] = atom
        ca = seen_names.get("CA")
        if ca is None:
            continue
        ca_coords.append([ca.pos.x, ca.pos.y, ca.pos.z])
        ca_resnums.append(res.seqid.num)
        for name, atom in seen_names.items():
            if atom.element.name != "H":
                heavy_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                heavy_resnums.append(res.seqid.num)
    if len(ca_coords) < 20:
        return None
    return {
        "chain": ch.name,
        "ca_coords": np.array(ca_coords),
        "ca_resnums": np.array(ca_resnums),
        "heavy_coords": np.array(heavy_coords),
        "heavy_resnums": np.array(heavy_resnums),
    }


# NOTE: an earlier draft of this script excluded crystallization-fusion
# inserts (T4-lysozyme, BRIL, etc.) via a resnum-floor heuristic (numbering
# >= 1000 => fusion partner). That heuristic was removed: it silently
# truncated genuine native residues for receptors whose own numbering
# legitimately exceeds 1000 with no gap at all (e.g. 7SF7 / ADGRL3, an
# adhesion GPCR with a long native N-terminal domain). Fusion partners are
# excluded implicitly and correctly by the TM-flag restriction in
# orient_structure()/main() below: they fall outside the UniProt-derived TM
# ranges, so they never enter the TM-restricted depth-binning profile,
# regardless of how they happen to be numbered. See VALIDATION.md.


# ---------------------------------------------------------------------------
# Orientation: three-tier procedure
# ---------------------------------------------------------------------------

def pca_normal(tm_coords):
    centroid = tm_coords.mean(axis=0)
    centered = tm_coords - centroid
    cov = np.cov(centered.T)
    evals, evecs = np.linalg.eigh(cov)
    return evecs[:, np.argmax(evals)], centroid


def orient_structure(acc, ca_coords, ca_resnums, pdb_id, n_terminus_side=None):
    """Return (depth_array_for_ca_coords, tm_flag, orient_method, orient_confidence).

    depth is signed distance (Angstrom) along the resolved membrane normal,
    with the convention that extracellular = positive, cytoplasmic =
    negative (verified where possible). tm_flag marks which residues fall
    within the UniProt-annotated transmembrane region -- this is used
    downstream to restrict the pivot search window, since structures with
    large extramembranous domains (e.g. Class C GPCR extracellular venus-
    flytrap domains, ion-channel cytoplasmic domains) otherwise produce a
    packing maximum far outside the membrane, which is never the intended
    "pivot".
    """
    tm_ranges, topo_domains = fetch_uniprot_tm_and_topology(acc)

    # Map residue numbers -> TM flag using UniProt TM ranges directly on
    # author-scheme resnums (works when CIF numbering matches UniProt
    # numbering; for offset numbering, align externally and supply an
    # already-TM-flagged coordinate set instead).
    tm_flag = np.array([any(s <= r <= e for s, e in tm_ranges) for r in ca_resnums])
    if tm_flag.sum() < 10:
        # Not enough resolved TM residues to get a stable normal; use all
        # CA atoms as a last resort for the PCA step only, and do not
        # restrict the pivot search window (no reliable TM boundary known).
        tm_flag_for_pca = np.ones(len(ca_resnums), dtype=bool)
        tm_flag_for_window = None
    else:
        tm_flag_for_pca = tm_flag
        tm_flag_for_window = tm_flag

    tm_coords = ca_coords[tm_flag_for_pca]
    normal, centroid = pca_normal(tm_coords)
    depth_all = (ca_coords - centroid) @ normal

    # Tier 1: topology verification
    ec_depths, ic_depths = [], []
    for i, r in enumerate(ca_resnums):
        for s, e, label in topo_domains:
            if s <= r <= e:
                if label == "Extracellular":
                    ec_depths.append(depth_all[i])
                elif label == "Cytoplasmic":
                    ic_depths.append(depth_all[i])
                break
    if len(ec_depths) >= 3 and len(ic_depths) >= 3:
        if np.mean(ec_depths) < np.mean(ic_depths):
            depth_all = -depth_all
        return depth_all, tm_flag_for_window, "topology_verified", "high"

    # Tier 2: PDBTM fallback
    tmatrix = fetch_pdbtm_matrix(pdb_id)
    if tmatrix is not None:
        depth_pdbtm = pdbtm_depth(ca_coords, tmatrix)
        return depth_pdbtm, tm_flag_for_window, "PDBTM_TMDET_verified", "high_PDBTM"

    # Tier 3: N/C-terminus fallback (requires caller-supplied side; we do
    # not guess it silently)
    if n_terminus_side in ("cytoplasmic", "extracellular"):
        n_term_depth = depth_all[0]
        want_negative = (n_terminus_side == "cytoplasmic")
        if (n_term_depth > 0) == want_negative:
            depth_all = -depth_all
        return depth_all, tm_flag_for_window, "Nterm_cyto_fallback", "low_Nterm_Cterm_disagree"

    # Unresolved: caller should exclude this structure from the primary set
    return depth_all, tm_flag_for_window, None, "unresolved"


# ---------------------------------------------------------------------------
# Pivot computation (unchanged from prior version)
# ---------------------------------------------------------------------------

def coordination_number_packing(coords_full, radius=10.0):
    """Per-atom coordination number (count of neighbors within `radius`,
    strict '<' cutoff), computed from the FULL set of atoms passed in
    (all Cα atoms of the chain, BEFORE restricting to the TM region) so
    that no genuine spatial neighbor is missed for residues near the
    TM-region boundary. This
    matches the manuscript's primary analysis exactly (see VALIDATION.md):
    the neighbor-counting pool is the whole chain; only the later
    depth-binning step (bin_and_find_pivot) restricts which residues
    contribute to -- and which bins are eligible for -- the reported pivot.
    """
    tree = cKDTree(coords_full)
    neighbor_counts = tree.query_ball_point(coords_full, r=radius, return_length=True) - 1
    return neighbor_counts.astype(float)


def all_atom_wcn(heavy_coords, capture_radius=15.0):
    """Per-heavy-atom Weighted Contact Number: WCN_i = sum_j 1/d_ij^2 over
    all other heavy atoms j within `capture_radius` (contributions beyond
    this are negligible under 1/d^2 weighting). Matches the manuscript's
    primary analysis exactly: computed over ALL heavy atoms of the chain
    (not just the TM region), consistent with coordination_number_packing
    above using the full chain as its neighbor pool.
    """
    tree = cKDTree(heavy_coords)
    pairs = tree.query_pairs(r=capture_radius, output_type="ndarray")
    n = len(heavy_coords)
    wcn = np.zeros(n)
    if len(pairs) > 0:
        d = np.linalg.norm(heavy_coords[pairs[:, 0]] - heavy_coords[pairs[:, 1]], axis=1)
        d = np.maximum(d, 0.5)  # avoid blow-up for near-overlapping/duplicate atoms
        contrib = 1.0 / (d ** 2)
        np.add.at(wcn, pairs[:, 0], contrib)
        np.add.at(wcn, pairs[:, 1], contrib)
    return wcn


def bin_and_find_pivot(depth_tm, values_tm, bin_width=3.0, readout_range=21.0):
    """Bin per-residue `values_tm` (already restricted to TM-flagged
    residues) by `depth_tm` into fixed-width bins over a FIXED range
    [-readout_range, +readout_range] (matching the manuscript's primary
    analysis exactly -- not a data-driven min/max range), and return
    (peak_depth, clarity). Pivot = the bin of maximum mean packing value
    (most tightly packed / most sterically constrained region). Bins with
    fewer than 3 contributing residues are excluded from the profile.
    clarity = peak value / minimum non-empty bin value (a peakedness
    statistic; NOT normalized by the profile median).
    """
    bins = np.arange(-readout_range, readout_range + 0.1, bin_width)
    mids, mm = [], []
    for i in range(len(bins) - 1):
        sel = (depth_tm >= bins[i]) & (depth_tm < bins[i + 1])
        if sel.sum() >= 3:
            mids.append((bins[i] + bins[i + 1]) / 2)
            mm.append(values_tm[sel].mean())
    if not mids:
        return None, None
    mids = np.array(mids)
    mm = np.array(mm)
    peak_depth = float(mids[np.argmax(mm)])
    clarity = float(mm.max() / mm.min()) if mm.min() > 0 else None
    return peak_depth, clarity


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/structure_manifest.csv")
    ap.add_argument("--cif-dir", default="data/cif")
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--out", default="results/pivot_depths.csv")
    args = ap.parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)
    packing_radius = cfg.get("coordination_number_packing", {}).get("radius_angstrom", 10.0)
    wcn_capture_radius = cfg.get("all_atom_wcn", {}).get("capture_radius_angstrom", 15.0)
    bin_width = cfg.get("depth_binning", {}).get("bin_width_angstrom", 3.0)
    readout_range = cfg.get("depth_binning", {}).get("readout_range_angstrom", 21.0)

    manifest = pd.read_csv(args.manifest)
    rows = []
    for _, row in manifest.iterrows():
        pdb_id, acc = row["pdb"], row.get("acc")
        chain_id = row["chain_id"] if "chain_id" in row and pd.notna(row.get("chain_id")) else None
        n_term_side = row["n_terminus_side"] if "n_terminus_side" in row and pd.notna(row.get("n_terminus_side")) else None
        cif_path = os.path.join(args.cif_dir, f"{pdb_id}.cif")
        if not os.path.exists(cif_path):
            print(f"{pdb_id}: SKIP (CIF not found at {cif_path})")
            continue

        parsed = get_ca_and_heavy(cif_path, chain_id=chain_id)
        if parsed is None:
            print(f"{pdb_id}: SKIP (chain not found / too few CA atoms)")
            continue

        # NOTE: crystallization-fusion partners (T4-lysozyme, BRIL, etc.) are
        # NOT excluded by a resnum-floor heuristic here. That approach was
        # tried and found to silently truncate genuine native residues for
        # receptors whose OWN numbering legitimately exceeds 1000 (e.g. 7SF7 /
        # ADGRL3, an adhesion GPCR with a long native N-terminal domain and no
        # numbering gap at all near 1000 -- see VALIDATION.md). Fusion
        # partners are excluded implicitly and correctly instead: they fall
        # outside the UniProt-derived TM ranges used for tm_flag, so they
        # never enter the TM-restricted depth-binning profile below,
        # regardless of how they are numbered.
        ca_coords, ca_resnums = parsed["ca_coords"], parsed["ca_resnums"]
        heavy_coords, heavy_resnums = parsed["heavy_coords"], parsed["heavy_resnums"]

        depth, tm_flag, orient_method, orient_confidence = orient_structure(acc, ca_coords, ca_resnums, pdb_id, n_term_side)

        if orient_confidence == "unresolved":
            print(f"{pdb_id}: orientation UNRESOLVED (no topology, no PDBTM, no N/C-term side given) -- excluded")
            rows.append({"pdb": pdb_id, "family": row.get("family"), "chain": parsed["chain"],
                         "orient_method": orient_method, "orient_confidence": orient_confidence,
                         "pivot_packing": np.nan, "clarity_packing": np.nan,
                         "pivot_wcn": np.nan, "clarity_wcn": np.nan})
            continue

        if tm_flag is None or tm_flag.sum() < 15:
            print(f"{pdb_id}: SKIP (fewer than 15 TM-flagged residues -- cannot build a reliable profile)")
            rows.append({"pdb": pdb_id, "family": row.get("family"), "chain": parsed["chain"],
                         "orient_method": orient_method, "orient_confidence": orient_confidence,
                         "pivot_packing": np.nan, "clarity_packing": np.nan,
                         "pivot_wcn": np.nan, "clarity_wcn": np.nan})
            continue

        # Coordination-number packing: computed over ALL Cα atoms of the
        # chain (ca_coords), then restricted to TM-flagged residues only
        # for depth-binning/pivot search over a FIXED [-21,21] A range.
        # This exactly matches the manuscript's primary analysis (see
        # VALIDATION.md) -- restricting only the SEARCH WINDOW post-hoc
        # (an earlier version of this script) does not reproduce it.
        coord_num_full = coordination_number_packing(ca_coords, radius=packing_radius)
        depth_ca_tm = depth[tm_flag]
        coord_num_tm = coord_num_full[tm_flag]
        pivot_packing, clarity_packing = bin_and_find_pivot(depth_ca_tm, coord_num_tm, bin_width=bin_width, readout_range=readout_range)

        # All-atom WCN: computed over ALL heavy atoms of the chain
        # (heavy_coords), aggregated to per-residue mean, then restricted
        # to TM-flagged residues (matched by resnum) for the same
        # fixed-range depth-binning as above.
        wcn_full = all_atom_wcn(heavy_coords, capture_radius=wcn_capture_radius)
        wcn_by_res = pd.Series(wcn_full).groupby(heavy_resnums).mean()
        depth_by_resnum = dict(zip(ca_resnums, depth))
        dp_list, wc_list = [], []
        for r in ca_resnums[tm_flag]:
            if r in wcn_by_res.index and r in depth_by_resnum:
                dp_list.append(depth_by_resnum[r])
                wc_list.append(wcn_by_res.loc[r])
        if len(dp_list) < 15:
            pivot_wcn, clarity_wcn = np.nan, np.nan
        else:
            pivot_wcn, clarity_wcn = bin_and_find_pivot(np.array(dp_list), np.array(wc_list), bin_width=bin_width, readout_range=readout_range)

        pp_str = f"{pivot_packing:.1f}" if pivot_packing is not None else "NA"
        pw_str = f"{pivot_wcn:.1f}" if pivot_wcn is not None and not (isinstance(pivot_wcn, float) and np.isnan(pivot_wcn)) else "NA"
        print(f"{pdb_id}: orient={orient_method} ({orient_confidence})  pivot_packing={pp_str}  pivot_wcn={pw_str}")
        rows.append({"pdb": pdb_id, "family": row.get("family"), "chain": parsed["chain"],
                     "orient_method": orient_method, "orient_confidence": orient_confidence,
                     "pivot_packing": pivot_packing, "clarity_packing": clarity_packing,
                     "pivot_wcn": pivot_wcn, "clarity_wcn": clarity_wcn})

    out_df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out_df.to_csv(args.out, index=False)
    print(f"\nSaved {len(out_df)} rows to {args.out}")
    n_unresolved = (out_df["orient_confidence"] == "unresolved").sum()
    print(f"Unresolved orientation (excluded from primary analysis): {n_unresolved}")


if __name__ == "__main__":
    main()
