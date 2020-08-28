"""
#!/usr/bin/env python
#title           :pyscript.py
#description     :This will merge UCSC expression, pheno, survival info into a big matrix as well as cancer by cancer
#author          :Rongbin Zheng
#date            :2020-08-29
#usage           :python prepare_data.py
#python_version  :3.6 

download expression data, phenotype, survival information from UCSC xena
#wget https://gdc.xenahubs.net/download/GDC-PANCAN.htseq_fpkm-uq.tsv.gz
#wget https://gdc.xenahubs.net/download/gencode.v22.annotation.gene.probeMap
#wget https://gdc.xenahubs.net/download/GDC-PANCAN.basic_phenotype.tsv.gz
#wget https://gdc.xenahubs.net/download/GDC-PANCAN.survival.tsv.gz
"""
import os
import sys
import pandas as pd
import numpy as np

# =========== FPKM-UQ version ============

phenotype_path = 'GDC-PANCAN.basic_phenotype.tsv'
phenotype = pd.read_csv(phenotype_path, sep='\t')

survial = pd.read_csv(
    'GDC-PANCAN.survival.tsv', sep='\t')

tcga_exp = pd.read_csv('GDC-PANCAN.htseq_fpkm-uq.tsv',
                       sep='\t', index_col=0)
# annotate gene name
gene_ann = pd.read_csv(
    'gencode.v22.annotation.gene.probeMap', sep='\t')

tcga_exp['gene_name'] = gene_ann.reindex(
    index=tcga_exp.index.tolist())['gene_name'].tolist()

c = pd.Series(collections.Counter(tcga_exp['gene_name']))
# unique
tcga_exp_uniq = tcga_exp[tcga_exp.gene_name.isin(c[c == 1].index.tolist())]
tcga_exp_uniq.index = tcga_exp.gene_name.tolist()
tcga_exp_uniq = tcga_exp_uniq.drop(['gene_name'], axis=1)
# multi
tcga_exp_multi = tcga_exp[tcga_exp.gene_name.isin(c[c > 1].index.tolist())]
tcga_exp_multi = tcga_exp_multi.groupby('gene_name').max()
tcga_exp = pd.concat([tcga_exp_multi, tcga_exp_uniq])
## merge phenotype
tcga_exp = pd.merge(tcga_exp.T, phenotype, left_index = True, right_on = 'sample')
tcga_exp.index = tcga_exp['sample'].tolist()
tcga_exp = tcga_exp.drop(['sample'], axis = 1)

del tcga_exp_uniq
del tcga_exp_multi
tcga_exp['sampleType'] = tcga_exp['sample_type'].tolist()
tcga_exp['sample_type'] = ['Normal' if x ==
                           'Solid Tissue Normal' else 'Tumor' for x in tcga_exp['sampleType'].tolist()]

# merge survial
tcga_exp = pd.merge(tcga_exp, survial, left_index=True, right_on='sample')
tcga_exp = tcga_exp.drop(['sample'], axis=1)

## output
os.mkdir('db')
tcga_exp.to_parquet('db/TCGA_XENA_exp_mat_with_pheno_survival.parquet')

# output by cancer type
os.mkdir('db/cancer')
for cancer in tcga_exp.project_id.unique():
	tcga_exp[tcga_exp['project_id'] == cancer].to_csv('db/cancer/%s.csv.gz', compression = 'gzip')



