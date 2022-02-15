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
    parser.add_argument("adata_file", help="normalized anndata input file")
    parser.add_argument("output_h5ad", help="output_file_name")
    parser.add_argument("hvgs", help="overdispersed_genes")
    parser.add_argument("runname", help="cNMF run name")
    args = parser.parse_args()
    anndata=args.adata_file
    output_h5ad=args.output_h5ad
    run_name=args.runname

    adata = sc.read(adata_file)