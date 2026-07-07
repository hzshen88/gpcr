#!/usr/bin/env python3
"""
Step 1: Fetch CIF structure files for the structure manifest.

Reads data/structure_manifest.csv (pdb, acc, chain_id, family, receptor,
group, excluded_orientation_unresolved) and downloads each structure from
RCSB PDB into data/cif/<PDB_ID>.cif. The manifest documents all 97
structures considered; 9 are flagged excluded_orientation_unresolved=True
and are excluded from the primary (n=88) analysis in
02_orient_and_compute_pivot.py -- see manuscript Methods, "Structure set".

Usage:
    python 01_fetch_structures.py [--manifest PATH] [--outdir PATH]
"""
import argparse
import os
import time
import urllib.request
import urllib.error

import pandas as pd


def fetch_cif(pdb_id: str, outdir: str, retries: int = 3, pause: float = 0.5) -> bool:
    dest = os.path.join(outdir, f"{pdb_id}.cif")
    if os.path.exists(dest):
        return True
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(url, dest)
            return True
        except urllib.error.URLError:
            time.sleep(pause)
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="data/structure_manifest.csv")
    ap.add_argument("--outdir", default="data/cif")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    manifest = pd.read_csv(args.manifest)

    failed = []
    for pdb_id in manifest["pdb"]:
        ok = fetch_cif(pdb_id, args.outdir)
        print(f"{pdb_id}: {'OK' if ok else 'FAILED'}")
        if not ok:
            failed.append(pdb_id)

    if failed:
        print(f"\n{len(failed)} structures failed to download: {failed}")
        print("Re-run this script to retry, or fetch these manually from https://www.rcsb.org/")
    else:
        print(f"\nAll {len(manifest)} structures downloaded successfully.")


if __name__ == "__main__":
    main()
