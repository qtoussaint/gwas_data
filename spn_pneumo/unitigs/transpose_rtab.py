#!/usr/bin/env python3
"""Transpose unitig-caller's native rtab (unitigs in rows, samples in columns)
into a samples x unitig-IDs 0/1 matrix.

Inputs:
  unitigs/spn_unitigs.rtab
      header: Unitig_sequence<TAB>sample1<TAB>sample2...
      body:   <unitig_sequence><TAB>0/1<TAB>0/1...

Outputs:
  spn_unitigs.matrix.tsv
      header: sample<TAB>unitig_1<TAB>unitig_2...<TAB>unitig_N
      one row per sample: <sampleID><TAB>0/1...   (sample IDs normalised)
  spn_unitigs.unitig_ids.tsv
      unitig_id<TAB>unitig_sequence   (maps each column ID back to its sequence)
"""
import os
import sys
import numpy as np

RTAB = sys.argv[1] if len(sys.argv) > 1 else "unitigs/spn_unitigs.rtab"
OUT_MATRIX = sys.argv[2] if len(sys.argv) > 2 else "unitigs/spn_unitigs.matrix.tsv"
OUT_IDS = sys.argv[3] if len(sys.argv) > 3 else "unitigs/spn_unitigs.unitig_ids.tsv"


def normalise(name):
    """Bifrost colour name -> bare sample ID (e.g. 6925_1#49)."""
    base = os.path.basename(name.strip())
    for suffix in (".contigs_velvet.fa", ".fa", ".fasta", ".fna"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def main():
    cols = []          # list of uint8 arrays, one per unitig (length = n_samples)
    sample_names = None
    n_unitigs = 0

    with open(RTAB) as fh, open(OUT_IDS, "w") as ids_out:
        header = fh.readline().rstrip("\n").split("\t")
        sample_names = [normalise(c) for c in header[1:]]
        n_samples = len(sample_names)

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            seq = parts[0]
            vals = parts[1:]
            if len(vals) != n_samples:
                sys.exit(f"Row width {len(vals)} != {n_samples} samples for unitig {n_unitigs+1}")
            n_unitigs += 1
            uid = f"unitig_{n_unitigs}"
            ids_out.write(f"{uid}\t{seq}\n")
            cols.append(np.array(vals, dtype=np.uint8))

    # (n_unitigs, n_samples) -> transpose -> (n_samples, n_unitigs)
    mat = np.vstack(cols).T if cols else np.empty((n_samples, 0), dtype=np.uint8)
    sys.stderr.write(f"matrix shape: {mat.shape[0]} samples x {mat.shape[1]} unitigs\n")

    unitig_ids = [f"unitig_{i+1}" for i in range(n_unitigs)]
    with open(OUT_MATRIX, "w") as out:
        out.write("sample\t" + "\t".join(unitig_ids) + "\n")
        for i, s in enumerate(sample_names):
            out.write(s + "\t" + "\t".join(map(str, mat[i].tolist())) + "\n")

    sys.stderr.write(f"wrote {OUT_MATRIX} and {OUT_IDS}\n")


if __name__ == "__main__":
    main()
