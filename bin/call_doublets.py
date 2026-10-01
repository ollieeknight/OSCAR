#!/usr/bin/env python3

import argparse

import pandas as pd
import scanpy as sc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5", required=True, help="cellbender filtered h5")
    ap.add_argument("--min-counts", type=int, required=True)
    ap.add_argument("--doublet-rate", type=float, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    adata = sc.read_10x_h5(args.h5)
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.scrublet(adata, expected_doublet_rate=args.doublet_rate)

    # The automatic threshold assumes a bimodal simulated-doublet histogram. When it is
    # unimodal the threshold lands in the tail and calls almost nothing, so call the
    # expected rate by quantile instead.
    score = adata.obs["doublet_score"]
    threshold = adata.uns["scrublet"]["threshold"]
    if (score > threshold).mean() < args.doublet_rate / 2:
        threshold = score.quantile(1 - args.doublet_rate)

    out = pd.DataFrame({"doublet_score": score, "is_gex_doublet": score > threshold}, index=adata.obs_names)
    out.index.name = "barcode"
    out.to_csv(args.out)


if __name__ == "__main__":
    main()
