import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt


def visualize_clustering(mat, linkage, out_path):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 10))
    dendro = hierarchy.dendrogram(linkage, no_plot=False, ax=ax1)
    g = sns.heatmap(mat.iloc[dendro['leaves'],dendro['leaves']], cmap='Blues',
     square=True, xticklabels=False, yticklabels=False, cbar=False, ax=ax2)
    for l in ['left','right','top','bottom']:
        g.spines[l].set_visible(True)
        g.spines[l].set_color('k')
    
    plt.save(f"{out_path}.plot.png")
    fig.close()


def main(input_matrix, input_matrix_ids, meta_path, outpath):
    new_meta_path = os.path.join(outpath, "metadata.merged.tsv") 
    indivs = np.loadtxt(input_matrix_ids, skiprows=0, dtype=str)
    rel_mat = np.loadtxt(input_matrix)
    mat = pd.DataFrame(rel_mat, index=indivs, columns=indivs)
    mat[mat < 0.4] = 0
    mat[np.isnan(mat)] = 0
    linkage = hierarchy.linkage(mat, method='complete', metric='correlation')
    cl = hierarchy.fcluster(linkage, 0.1, criterion='distance')
    clusters = pd.DataFrame({'ds_number': mat.index, 'genotype_cluster': cl}).sort_values(
        by='genotype_cluster')
    
    metadata = pd.read_table(meta_path, header=0, dtype={'indiv_id': object})
    metadata = metadata.rename(columns={'indiv_id': 'ds_number'})
    metadata = metadata.merge(clusters, on='ds_number').sort_values(by='genotype_cluster')
    metadata = metadata.rename(columns={'genotype_cluster': 'indiv_id'})
    metadata['indiv_id'] = 'INDIV_' + metadata['indiv_id'].astype(str).str.zfill(4)
    metadata.to_csv(new_meta_path, header=True, index=False, sep='\t')
    visualizations_path = os.path.join(outpath, 'clustering.png')
    visualize_clustering(mat, linkage, out_path=visualizations_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count tags by allele")
    parser.add_argument("--matrix", type=str, help="Result of plink2 --king with suffix .king")

    parser.add_argument("--matrix-ids", type=str,
						help="Result of plink2 --king with suffix .king.id")

    parser.add_argument("--meta-file", type=str, 
						help="Path to meta file")

    parser.add_argument("--outpath", type=str, 
						help="Output path to directory with meta file and")

    args = parser.parse_args()

    main(args.matrix, args.matrix_ids, args.meta_file, args.outpath)