#!/usr/bin/env python3

def read_gmt(gmt_file):
    genesets = {}
    with open(gmt_file, 'r') as f:
        unsplit = f.read().split('\n')
        for sig in unsplit:
            gene_arr = sig.split('\t')
            if gene_arr[0] != '':
                # first entry is name, skip description, rest are genes
                genesets[gene_arr[0]] = gene_arr[2:]
    return genesets

def write_gmt(gmt, outfile):
    with open(outfile, 'w') as f:
        for k in gmt.keys():
            f.write(k + '\t_\t')
            f.write('\t'.join(gmt[k]) + '\n')

