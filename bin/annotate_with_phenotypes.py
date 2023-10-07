import os
import pandas as pd
import glob
import sys
from tqdm import tqdm
from functools import reduce


tqdm.pandas()
phenotype_db_names = ['grasp', 'ebi', 'clinvar', 'phewas', 'finemapping']


def pack(arr):
    return '\t'.join(map(str, arr)) + '\n'


def parse_grasp(filepath):
    print('Parsing Grasp')
    with open(filepath) as file:
        r = file.readlines()
    result = {'ID': [], 'grasp': []}
    for line in tqdm(list(r)):
        a = line.strip('\n').split('\t')
        result['ID'].append(f"rs{a[0]}")
        result['grasp'].append(a[1])
    return pd.DataFrame(result).groupby('ID').agg(
        grasp=('grasp', arr_to_str)
    ).reset_index()


def parse_finemapping(filepath):
    result = {'ID': [], 'finemapping': []}
    print('Parsing Finemapping')
    with open(filepath) as file:
        f = file.readlines()
    for k, line in enumerate(tqdm(list(f))):
        if k == 0:
            continue
        a = line.strip('\n').split(',')

        result['ID'].append(a[2])
        result['finemapping'].append(a[0])
    return pd.DataFrame(result).groupby('ID').agg(
        finemapping=('finemapping', arr_to_str)
    ).reset_index()

def parse_clinvar(filepath):
    result = {'ID': [], 'clinvar': []}
    print('Parsing ClinVar')
    with open(filepath, 'r') as file:
        f = file.readlines()
    for k, line in enumerate(tqdm(list(f))):
        if k == 0:
            continue
        a = line.strip('\n').split('\t')
        if 'pathogenic' not in a[6].lower() and 'risk factor' not in a[6].lower():
            continue

        for ph in a[13].split(';'):
            if ph in ('not provided', 'not specified'):
                continue
            result['ID'].append(f"rs{a[9]}")
            result['clinvar'].append(ph)
    return pd.DataFrame(result).groupby('ID').agg(
        clinvar=('clinvar', arr_to_str)
    ).reset_index()


def parse_ebi(filepath):
    result = {'ID': [], 'ebi': []}
    print('Parsing EBI')
    with open(filepath, 'r') as file:
        f = file.readlines()
    for k, line in enumerate(tqdm(list(f))):
        if k == 0:
            continue
        try:
            a = line.strip('\n').split('\t')
            result['ID'].append(f"rs{a[23]}")
            result['ebi'].append(a[7])
        except ValueError:
            continue
    return pd.DataFrame(result).groupby('ID').agg(
        ebi=('ebi', arr_to_str)
    ).reset_index()


def parse_phewas(filepath):
    result = {'ID': [], 'phewas': []}
    print('Parsing phewas')
    with open(filepath, 'r') as file:
        f = file.readlines()
    for k, line in enumerate(tqdm(list(f))):
        if k == 0:
            continue
        a = line[line.find('"rs'):].split('",')
        ph = a[1][1:]
        rs_id = f"rs{int(a[0][3:])}"
        result['ID'].append(rs_id)
        result['phewas'].append(ph)
    return pd.DataFrame(result).groupby('ID').agg(
        phewas=('phewas', arr_to_str)
    ).reset_index()


def parse_gtex(snps_pos, qtlfiles, transqtl):

    cis = {'#chr': [], 'end': [], 'cis': []}
    trans = {'#chr': [], 'end': [], 'trans': []}
    print('Number of cis-eQTL files:', len(qtlfiles))

    for qtlfile in tqdm(list(qtlfiles)):
        with open(qtlfile) as qfile:

            tis = qtlfile[qtlfile.rfind('/') + 1:qtlfile.find('.v8.')]
            for line in qfile:

                if line.startswith('variant_id'):
                    tit = line.strip('\n').split('\t')
                    titlen = len(tit)
                    continue

                a = line.strip('\n').split('\t')
                a = {tit[x]: a[x] for x in range(titlen)}

                chrom, end = a['variant_id'].split('_')[:2]
                cis['#chr'].append(chrom)
                cis['end'].append(int(end))
                cis['cis'].append(a['gene_id'])

    print('Starting to parse transqtl')
    with open(transqtl) as trfile:
        for line in trfile:
            if line.startswith('tissue_id'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}

            chrom, end = a['variant_id'].split('_')[:2]
            tis = a['tissue_id']
            gen = a['gene_id']
            trans['#chr'].append(chrom)
            trans['end'].append(int(end))
            trans['trans'].append(gen)
    trans = pd.DataFrame(trans).groupby(['#chr', 'end']).agg(
        trans=('trans', arr_to_str)
    ).reset_index()
    cis = pd.DataFrame(cis).groupby(['#chr', 'end']).agg(
        cis=('cis', arr_to_str)
    ).reset_index()
    return snps_pos.merge(cis, how='left').merge(trans, how='left')


def arr_to_str(arr):
    non_nans = [remove_phen_name_punctuation(x) for x in arr if x is not None]
    if len(non_nans) == 0:
        return ""
    return '|'.join(sorted(non_nans))


def remove_phen_name_punctuation(phenotype_name):
    return phenotype_name.lower().replace("'", '').replace('_', ' ')

def parse_dbs(snps_positions, grasp, ebi, clinvar, fm, phewas):
    g = parse_grasp(grasp)
    e = parse_ebi(ebi)
    c = parse_clinvar(clinvar)
    p = parse_phewas(phewas)
    f = parse_finemapping(fm)
    dfs = [snps_positions, g, e, c, p, f]
    return reduce(lambda left, right: pd.merge(
        left, right,
        on=['ID'],
        how='left'
    ),
    dfs
    )

def main(phenotypes_dir, snps_path, out_path):
    print('Reading files')
    snps_positions = pd.read_table(snps_path)
    grasp = os.path.join(phenotypes_dir, 'pheno', 'grasp_pheno.tsv')
    ebi = os.path.join(phenotypes_dir, 'pheno', 'gwas_catalog.tsv')
    clinvar = os.path.join(phenotypes_dir, 'pheno', 'variant_summary.txt')
    fm = os.path.join(phenotypes_dir, 'pheno', 'finemapping.csv')
    phewas = os.path.join(phenotypes_dir, 'pheno', 'phewas-catalog.csv')

    qtlfiles = glob.glob(os.path.join(phenotypes_dir, 'eqtl', 'signif', '*.txt'))
    transqtl = os.path.join(phenotypes_dir, 'eqtl', 'GTEx_Analysis_v8_trans_eGenes_fdr05.txt')
    print('Started parsing DBs')
    phen_df = parse_dbs(snps_positions, grasp, ebi, clinvar, fm, phewas)

    print('Started parsing GTEX')
    phen_df = parse_gtex(phen_df, qtlfiles, transqtl)
    print('Parsing finished')
    phen_df.to_csv(out_path, sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])