""" Parse results from ps_scan 
    OUTPUT:
    protein start stop domains seq ProfileScan
"""
import sys

with open(sys.argv[1]) as f:
    for line in f:
        sp = line.strip().split('\t')
        print(sp[0] + '\t'
              + sp[1] + '\t' + sp[2] + '\t'
              + sp[3] + '\tseq\tProfileScan')
