import sys, os

import pandas as pd
import numpy as np
from umi_tools import UMIClusterer

from tqdm import tqdm


if __name__ == "__main__":
    # input is barcode_table.txt from proc_barcodes.sh
    infile = sys.argv[1]
    outfile = os.path.splitext(infile)[0] + '.EC.txt'
    
    read_table = pd.read_table(infile)
    
    print('correcting cell and sgRNA UMIs')
    clusterer = UMIClusterer(cluster_method='directional')
    cstack = []
    gstack = []
    for k, r in tqdm(read_table.groupby(['cell', 'sgRNA'])):
        if len(r) == 1:
            cstack.append(r['cell_UMI'])
            gstack.append(r['sgUMI'])
        else:
            # cell UMI
            umi = r['cell_UMI'].value_counts()
            umi = {str.encode(k): v for k, v in zip(umi.index, umi.values)}
            clusters = clusterer(umi, threshold=1)
            # first UMI in cluster is always the most abundant one
            corrected = r['cell_UMI'].apply(
                lambda x: clusters[np.where([str.encode(x) in _ for _ in clusters])[0][0]][0].decode())
            cstack.append(corrected)
    
            # sgUMI
            umi2 = r['sgUMI'].value_counts()
            umi2 = {str.encode(k): v for k, v in zip(umi2.index, umi2.values)}
            clusters = clusterer(umi2, threshold=1)
            corrected2 = r['sgUMI'].apply(
                lambda x: clusters[np.where([str.encode(x) in _ for _ in clusters])[0][0]][0].decode())
            gstack.append(corrected2)
    cstack = pd.concat(cstack)
    gstack = pd.concat(gstack)
    read_table['cell_UMI_EC'] = cstack.loc[read_table.index]
    read_table['sgUMI_EC'] = gstack.loc[read_table.index]
    
    n_corrected1 = (read_table['cell_UMI'] != read_table['cell_UMI_EC']).sum()
    n_corrected2 = (read_table['sgUMI'] != read_table['sgUMI_EC']).sum()
    print('{:.4f}% ({}/{}) cell UMIs corrected'.format(n_corrected1 / len(read_table) * 100, n_corrected1, len(read_table)))
    print('{:.4f}% ({}/{}) sgRNA UMIs corrected'.format(n_corrected2 / len(read_table) * 100, n_corrected2, len(read_table)))
    
    read_table.to_csv(outfile, sep='\t', index=False)

