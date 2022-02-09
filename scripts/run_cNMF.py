import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scanpy as sc

import argparse

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="""
    This script will run cNMF on counts data cellXgene
    matrix with cells as rows and genes as columns""")
    parser.add_argument("adata_file", help="anndata input file")
    parser.add_argument("outdir", help="output directory")
    parser.add_argument("output_h5ad", help="output_file_name")
    parser.add_argument("runname", help="cNMF run name")
    args = parser.parse_args()
    anndata=args.adata_file
    OUTDIR=args.outdir
    output_h5ad=args.output_h5ad
    run_name=args.runname

    adata = sc.read(anndata)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=200)
    sc.pp.filter_genes(adata, min_cells=3) #

    numiter=20
    numhvgenes=2000
    output_directory = OUTDIR

    os.chdir( OUTDIR )
    sc.write(output_h5ad,adata)
