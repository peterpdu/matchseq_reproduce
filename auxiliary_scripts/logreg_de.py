#!/usr/bin/env python3.6

"""
Logistic Regression differential gene expression
v1.0

Purpose:    does a LRT per-gene logistic regression to classify cell labels to find DEGs between two groups
Input:      scanpy h5ad file with cell identity labels; comparisons to perform
Output:     table(s) of differential genes with various stats

Version history:
v1.0: from Ntranos et al., 2019 [Nat Methods] from Pachter Lab

Last modified by Peter Du on Feb 2020
"""


import sys, os, argparse
import numpy as np
import pandas as pd

import scanpy as sc
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
from scipy import stats
from scipy.sparse import issparse
from statsmodels.stats.multitest import multipletests

from tqdm import tqdm
from datetime import datetime


def getParser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     add_help=True)
    required = parser.add_argument_group('required')
    opt = parser.add_argument_group('optional')
    required.add_argument('--h5ad',
                          type=str,
                          help='scanpy h5ad file',
                          required=True)
    required.add_argument('--compare',
                          type=str,
                          help='groups to compare; can be two column file or comma-separated names '
                               '(e.g. CD8_1,CD82 CD4,Treg ...); if only one value supplied, assumed one-vs-rest',
                          required=True,
                          nargs='+')

    opt.add_argument('--label',
                     type=str,
                     help='name of comparison label in adata object; default=celltype',
                     default='celltype')
    opt.add_argument('--group',
                     type=str,
                     help='group by this label and perform DE within these groupings; default=None',
                     default=None)
    opt.add_argument('--layer',
                     type=str,
                     help='which layer to use; default=adata.X',
                     default=None)
    return parser



def logreg_de(adata, comparison, labelname):
    """
    Parameters
    ----------
    adata: AnnData object
    comparison: one or two length iterable; if one length, one vs else; if two length, group A vs B
    labelname: adata.obs column to use for comparison label

    Returns: differential expression results table
    -------

    """
    if len(comparison) == 2:
        x1 = adata[adata.obs[labelname] == comparison[0]].X
        x2 = adata[adata.obs[labelname] == comparison[1]].X
    else:
        x1 = adata[adata.obs[labelname] == comparison[0]].X
        x2 = adata[adata.obs[labelname] != comparison[0]].X

    if issparse(adata.X):
        x1 = np.array(x1.todense())
        x2 = np.array(x2.todense())
    n1 = x1.shape[0]
    n2 = x2.shape[0]
    if (n1 == 0) or (n2 == 0):
        print('Error: no cells found for a group')
        return None
    if ((n1 / n2) > 5) or (n1 / n2) < 0.2:
        print('Warning: greater than 5-fold difference in cell counts between groups')

    local_labels = np.concatenate((np.ones(n1), np.zeros(n2)), axis=0)
    logr = LogisticRegression(solver='lbfgs')
    p_of_1 = n1 / float(n1 + n2)
    llnull = (n1 + n2) * (p_of_1 * np.log(p_of_1) + (1 - p_of_1) * np.log(1 - p_of_1))

    results = {'pvals': [], 'log_loss': [],
               'll_full': [], 'll_reduced': [],
               'mean': []}
    for i, gn in tqdm(enumerate(adata.var.index), total=len(adata.var)):
        local_x1 = x1[:, i]
        local_x2 = x2[:, i]
        mean = local_x1.mean() - local_x2.mean()
        cat = np.concatenate([local_x1, local_x2]).reshape(-1, 1)
        # c = c - np.mean(c, axis=0)
        try:
            logr.fit(cat, local_labels)
        except:
            continue
        pred = np.array(logr.predict_proba(cat)[:, 1])
        gene_score = log_loss(local_labels, pred)
        llf = -gene_score * (n1 + n2)
        llr = llf - llnull
        llr_pval = stats.chi2.sf(2 * llr, 1)

        results['pvals'].append(llr_pval)
        results['log_loss'].append(gene_score)
        results['ll_full'].append(llf)
        results['ll_reduced'].append(llr)
        results['mean'].append(mean)

    outdf = pd.DataFrame.from_dict(results)
    outdf.index = adata.var.index

    # Benjamini-Hochberg FDR correction
    outdf['fdr'] = multipletests(outdf['pvals'], method='fdr_bh')[1]

    outdf['signed_logfdr'] = np.sign(outdf['mean']) * -1 * np.log10(outdf['fdr'])
    outdf = outdf.sort_values('log_loss', ascending=True)

    return outdf


def main():
    args = getParser().parse_args()
    timestamp = str(datetime.now().strftime("%y-%m-%d-%H-%M-%S"))
    print('Reading h5ad')
    adata = sc.read_h5ad(args.h5ad)

    outdir = os.path.join(os.path.dirname(args.h5ad), 'logreg_de_{}'.format(timestamp))
    os.mkdir(outdir)

    labelname = args.label
    groupname = args.group

    layer_use = args.layer
    if layer_use is not None:
        print(f'using layer: {layer_use}')
        adata.X = adata.layers[layer_use]
    else:
        print('using adata.X')

    if groupname is not None:
        assert groupname in adata.obs.columns
        adata.obs[groupname] = adata.obs[groupname].astype(str)
        adata.obs[groupname] = [x.replace('/', '_') for x in adata.obs[groupname]]

    assert labelname in adata.obs.columns
    adata.obs[labelname] = adata.obs[labelname].astype(str)
    adata.obs[labelname] = [x.replace('/', '_') for x in adata.obs[labelname]]

    valid_comp_groups = adata.obs[labelname].unique()
    print('\nComparison groups found')
    for g in valid_comp_groups:
        print(g)
    print()

    if len(args.compare) == 1 and os.path.exists(args.compare[0]):
        print('Comparison file input detected')
        comparisons = pd.read_csv(args.compare[0], sep='\t').values
    else:
        comparisons = []
        for entry in args.compare:
            entry = entry.split(',')
            # check validity
            if not all([x in valid_comp_groups for x in entry]):
                print('{} contains invalid comparison group(s)'.format(entry))
                continue
            comparisons.append(entry)

    if len(comparisons) == 0:
        print('No valid comparisons found')
        sys.exit(1)

    # do DE by groups
    if groupname is not None:
        groups = adata.obs[groupname].unique()
        for g in groups:
            print('Group {}'.format(g))
            subset = adata[adata.obs[groupname] == g].copy()
            for c in comparisons:
                print('Analyzing {}'.format(c))
                outdf = logreg_de(subset, c, labelname)
                if outdf is not None:
                    outpath = os.path.join(outdir, '_'.join([g] + c + ['lrLRT.txt']))
                    print('Writing results to {}'.format(outpath))
                    outdf.to_csv(outpath, sep='\t')
            print()

    # do DE
    for c in comparisons:
        print('Analyzing {}'.format(c))
        outdf = logreg_de(adata, c, labelname)
        if outdf is not None:
            outpath = os.path.join(outdir, '_'.join(c + ['lrLRT.txt']))
            print('Writing results to {}'.format(outpath))
            outdf.to_csv(outpath, sep='\t')


if __name__ == '__main__':
    main()
