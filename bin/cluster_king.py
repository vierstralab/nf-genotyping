import argparse
from curses import meta
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
    
    plt.savefig(out_path)
    plt.close(fig)


def main(mat, genotyping_meta, outdir, min_hets=100):
    good_ids = genotyping_meta[genotyping_meta['nHets'] > min_hets].index.tolist()
    mat = mat.loc[good_ids, good_ids]

    linkage = hierarchy.linkage(mat, method='complete', metric='correlation')
    cl = hierarchy.fcluster(linkage, 0.1, criterion='distance')
    clusters = pd.DataFrame({'indiv_id': mat.index, 'genotype_cluster': cl}).sort_values(
        by='genotype_cluster')
    clusters['new_id'] = 'INDIV_' + metadata['genotype_cluster'].astype(str).str.zfill(5)

    
    metadata = metadata.merge(clusters,
        on='indiv_id', how='outer').sort_values(by='genotype_cluster')
    metadata.rename(columns={'indiv_id': 'old_indiv_id'}, inplace=True)
    metadata.rename(columns={'new_id': 'indiv_id'}, inplace=True)
    metadata['indiv_id'] = 'INDIV_' + metadata['indiv_id'].astype(str).str.zfill(5)
    metadata.drop(columns=['genotype_cluster'], inplace=True)
    
    metadata.to_csv(
        os.path.join(outdir, "metadata.clustered.tsv"),
        index=False, sep='\t'
    )

    visualize_clustering(
        mat, linkage,
        os.path.join(outdir, 'clustering.png')
    )


def read_genotype_stats(bcftools_stats):
    stats = pd.read_table(bcftools_stats)
    stats = stats.iloc[:,2:]
    stats.columns = stats.columns.str.replace(r"(\[.*\])","")
    stats.rename(columns={'sample': 'indiv_id'}, inplace=True)
    stats.set_index('indiv_id', inplace=True)
    return stats


def read_plink(plink_path):
    indivs = np.loadtxt(f'{plink_path}.king.id', skiprows=0, dtype=str)
    rel_mat = np.loadtxt(f'{plink_path}.king')
    mat = pd.DataFrame(rel_mat, index=indivs, columns=indivs)
    mat[mat < 0.4] = 0
    mat[np.isnan(mat)] = 0
    return mat

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count tags by allele")
    parser.add_argument('plink', type=str, help="Path to plink2 --king result files without (.king/.king.id) suffix")
    parser.add_argument('metadata', type=str, 
						help="Path to metadata file")
    
    parser.add_argument('stats', type=str, help='Path to bcftools PSC stats file')

    parser.add_argument('--outpath', type=str, 
						help="Path to directory to save updated metafile and visualizations", default='./')
    
    parser.add_argument('--min-hets',
						help="Minimal number of heterozygotes per sample to cluster", type=int, default=100)
    

    args = parser.parse_args()
    matrix = read_plink(args.plink)
    stats = read_genotype_stats(args.stats)
    metadata = pd.read_table(args.metadata, 
        dtype={'indiv_id': str}).join(stats)
    main(matrix, metadata, args.outpath, min_hets=args.min_hets)
    