import argparse
import pandas as pd
import matplotlib.pyplot as plt

from genome_tools.data.genotypes import read_king_matrix, read_bcftools_psc, cluster_relatedness
from genome_tools.plotting.relatedness import relatedness_heatmap


def main(mat, stats, genotyping_meta, min_variants=1000):
    good_ids = stats.query(f'(nHets + nNonRefHom) >= {min_variants}')['indiv_id'].tolist()
    mat = mat.loc[good_ids, good_ids]
    clusters, linkage = cluster_relatedness(mat)

    result = genotyping_meta.merge(
        clusters[['indiv_id', 'clustered_indiv_id']],
        how='left'
    ).rename(
        columns={'indiv_id': 'old_indiv_id', 'clustered_indiv_id': 'indiv_id'},
    )
    return mat, linkage, result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Count tags by allele")
    parser.add_argument('plink', type=str, help="Path to plink2 --king result files without (.king/.king.id) suffix")
    parser.add_argument('metadata', type=str, 
						help="Path to metadata file")
    
    parser.add_argument('stats', type=str, help='Path to bcftools PSC stats file')

    parser.add_argument('--outpath', type=str, 
						help="Path to directory to save updated metafile and visualizations", default='./')
    
    parser.add_argument('--min-variants',
        help="Minimal number of variants (Het + NonRefHom) per sample to include in clustering",
        type=int, default=1000
    )

    args = parser.parse_args()
    matrix = read_king_matrix(args.plink)
    stats = read_bcftools_psc(args.stats)
    
    metadata = pd.read_table(args.metadata, dtype={'indiv_id': str})
    
    mat, linkage, gen_meta = main(matrix, stats, metadata, min_variants=args.min_variants)

    gen_meta.to_csv(
        f"{args.outpath}/metadata.clustered.tsv",
        index=False,
        sep='\t'
    )

    fig, (ax_dendrogram, ax) = plt.subplots(2, 1, figsize=(6, 10))
    relatedness_heatmap(mat, linkage, ax=ax, ax_dendrogram=ax_dendrogram)
    fig.savefig(f"{args.outpath}/clustering.pdf", bbox_inches='tight', transparent=True)
    plt.close(fig)

    