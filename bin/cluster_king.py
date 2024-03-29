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
    
    plt.savefig(out_path)
    plt.close(fig)

def get_clusters(mat):
    linkage = hierarchy.linkage(mat, method='complete', metric='correlation')
    cl = hierarchy.fcluster(linkage, 0.1, criterion='distance')

    clusters = pd.DataFrame({'indiv_id': mat.index, 'cluster_id': cl})

    frequency = clusters['cluster_id'].value_counts().reset_index().reset_index(names='genotype_cluster')
    clusters = clusters.merge(frequency, on='cluster_id').sort_values(by='count', ascending=False)

    clusters['new_id'] = 'INDIV_' + (clusters['genotype_cluster'] + 1).astype(str).str.zfill(4)
    return clusters, linkage, cl


def main(mat, stats, genotyping_meta, outdir, min_hets=100):
    good_ids = stats.query(f'nHets >= {min_hets}')['indiv_id'].tolist()
    clusters, _, _ = get_clusters(mat.loc[good_ids, good_ids])

    genotyping_meta = genotyping_meta.merge(clusters[['indiv_id', 'new_id']], how='left')

    genotyping_meta.rename(
        columns={'indiv_id': 'old_indiv_id', 'new_id': 'indiv_id'},
        inplace=True
    )
    
    genotyping_meta.to_csv(
        os.path.join(outdir, "metadata.clustered.tsv"),
        index=False, sep='\t'
    )

    visualize_clustering(
        mat, linkage,
        os.path.join(outdir, "clustering.png")
    )


def read_genotype_stats(bcftools_stats):
    stats = pd.read_table(bcftools_stats)
    stats = stats.iloc[:,2:]
    stats.columns = stats.columns.str.replace(r"(\[.*\])","", regex=True)
    stats.rename(columns={'sample': 'indiv_id'}, inplace=True)
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
    
    metadata = pd.read_table(args.metadata, dtype={'indiv_id': str})
    
    main(matrix, stats, metadata, args.outpath, min_hets=args.min_hets)
    