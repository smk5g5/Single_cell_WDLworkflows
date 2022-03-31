import os
import pandas as pd
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
import scanpy as sc

import argparse

#This script will be updated to work with WDL

if __name__ == "__main__":
	parser=argparse.ArgumentParser(description="""
	This script will run cNMF on counts data cellXgene
    matrix with cells as rows and genes as columns""")

	parser.add_argument("Counts", help="Count input file")
	parser.add_argument("output_h5ad", help="output_file_name")

	args = parser.parse_args()

	COUNTS=args.Counts
	output_h5ad=args.output_h5ad
	f_exprMat = str(COUNTS)
	print("f_exprMat ",f_exprMat)
	adata = sc.read_text( f_exprMat, delimiter='\t', first_column_names=True )
	adata = adata.transpose()
	sc.write(output_h5ad, adata)

	