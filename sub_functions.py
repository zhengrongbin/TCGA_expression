"""
#!/usr/bin/env python
#author          :Rongbin Zheng
#date            :2020-08-29
#python_version  :3.6 
"""
import basic_tumor_normal as basic
import scipy.cluster.hierarchy as sch
from scipy import linalg
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter, CoxPHFitter
import os
import sys
import math
import re
import traceback
import pandas as pd
import numpy as np
import collections
import pickle as pk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy.stats.stats import pearsonr, spearmanr
from scipy.stats import wilcoxon, ranksums
from scipy.stats import stats
import argparse
from argparse import RawDescriptionHelpFormatter
import corr_gene_tumor as corr_gene
import advanced


def info(mesg):
    # output message
    os.system('echo "++++ %s"' % mesg)


def _read_in_file(path, idc):
    """
    read in expression matrix by files
    """
    info('read in file %s' % path)

    if not os.path.exists(path):
        info('file path not exist: %s' % path)
        sys.exit(1)
    try:
        if path.endswith('csv.gz'):
            mat = pd.read_csv(path, compression='gzip', index_col=0)
        elif path.endswith('.parquet'):
            mat = pd.read_parquet(path)
        else:
            mat = pd.read_csv(path, sep='\t', index_col=0)
    except:
        traceback.print_exc(file=sys.stderr)  # maybe the file type problem
        sys.exit(1)
    # TARGET-RT, too few sample is avaliable
    mat = mat[~mat.project_id.isin(['TARGET-RT'])]
    # check file title
    if 'project_id' not in mat.columns.tolist():
        info('project_id not in column names')
        sys.exit(1)
    if 'sample_type' not in mat.columns.tolist():
        info('sample_type is not in columns')
        sys.exit(1)
    # specify to needed genes:
    # the gene not in matrix columns
    diffgene = list(set(idc) - set(mat.columns.tolist()))
    if diffgene:
        info('these genes %s are not in the expression matrix of this cancer, skip %s' % (
            str(diffgene), str(path)))
        # return(pd.DataFrame())  # return a empty dataframe
    return (mat)


def _check_cancer(cancer, all_exist, dpath):
    # get cancer list that expression matrix is in folder
    if cancer == 'all':
        info('will read in all cancer types in %s' % dpath)
        return(all_exist)
    else:
        cancer = cancer.split(',') if ',' in cancer else [cancer]
        need_cancer = list(set(cancer) & set(all_exist))
        diff_cancer = list(set(cancer) - set(all_exist))
        if not need_cancer:  # no cancer type matched
            info(
                '-c is not all and none of given {0} in existed data {1}'.format(str(cancer), str(all_exist)))
            sys.exit(1)
        if diff_cancer:
            info('some cancer type %s is not existed in the folder %s' %
                 (str(diff_cancer, dpath)))
        return(need_cancer)


def _read_in_folder(dpath, cancer, idc):
    """
    read in user specified cancer exp mat to save memory
    """
    info('input folder is given, read in by cancer types one by one')

    all_exist = {x.split('.')[0]: os.path.join(dpath, x)
                 for x in os.listdir(dpath)}
    need_cancer = _check_cancer(cancer, list(all_exist.keys()), dpath)
    mat = pd.DataFrame()
    for cancer in need_cancer:
        info('read in exp mat of %s' % cancer)
        # print(idc)
        mat_tmp = _read_in_file(all_exist[cancer], idc)
        mat = pd.concat([mat, mat_tmp])  # rbind
        del mat_tmp
    return(mat)


def _need_genes(config):
    """extract needed genes, so only load data of these genes"""
    need_genes = []
    for t in ['gene', 'gene1', 'gene2']:
        if (t in config.keys()) and config[t]:
            need_genes.append(config[t])
    if ('adj_gene' in config.keys()) and config['adj_gene']:
        if config['adj_gene'] == 'CTL':
            need_genes.extend(['CD8A', 'CD8B', 'PRF1', 'GZMA', 'GZMB'])
        else:
            need_genes.append(config['adj_gene'])
    if ('protein_gene' in config.keys()) and config['protein_gene']:
        need_genes.extend(config['protein_gene'])
    return(need_genes)


def _read_in(config):
    """
    read in expression file by file or folder
    """
    # specify needed genes
    need_genes = _need_genes(config)
    idc = need_genes+['sample', 'project_id', 'sample_type', 'sampleType',
                      'OS', '_PATIENT', 'OS.time']

    if config['fpath']:
        # user gives file path
        mat = _read_in_file(config['fpath'], idc)
    elif config['dpath']:
        # user gives foder path where file was saved by cancer type
        mat = _read_in_folder(config['dpath'], config['cancer'], idc)
    else:
        info('Please set -i or -d')
        sys.exit(1)
    info('read in exp successfully')
    # check mat
    if mat.shape[0] == 0:
        info('No expression data loaded, please check reference files and given gene names')
        sys.exit(1)
    # check CTL
    if 'adj_gene' in config.keys() and config['adj_gene'] == 'CTL':
        mat['CTL'] = mat[['CD8A', 'CD8B', 'GZMB', 'GZMA', 'PRF1']].T.mean()
    return(mat)


def _basic(args, config):
    """basic job linker"""
    info('running on basic function')
    config['distribution'], config['boxplot'], config['survival'] = args.distribution, args.boxplot, args.survival
    config['metastatic'] = args.metastatic
    config['gene'] = args.gene_name
    if (not config['gene']) or (not config['distribution'] and not config['boxplot'] and not config['survival']) and not config['metastatic']:
        info('gene name or subplot argument not be specified')
        sys.exit(1)

    exp_mat = _read_in(config)
    config['exp_mat'] = exp_mat
    if config['distribution']:
        basic._distribution(config)
     # boxplot
    if config['boxplot']:
        basic._tumor_normal(config)

    # survial analysis
    if config['survival']:
        basic._run_survival(config)

    if config['metastatic']:
        basic._metatstic_boxplot(config)


def _corr(args, config):
    """correlation of given genes"""
    info('running on correlation function')
    config['gene1'], config['gene2'] = args.gene1, args.gene2
    config['adj_gene'] = args.adjust
    if config['gene1'] and config['gene2']:
        exp_mat = _read_in(config)
        config['exp_mat'] = exp_mat
    else:
        info('gene1 or gene2 error')
        sys.exit(1)
    corr_gene._corr_tumor(config)


def _gene_ann(gene_ann_path):
    """
    read in gene annotation, and get protein coding gene list
    since _advance() is design for calculating correlation of given gene against for all protein coding genes,
    so just need protein coding genes and or the given non-coding gene
    gene annotation could download from GENCODE, and clear into matrix where "gene_name" and "gene_type" are in the columns
    """
    gene_ann = pd.read_csv(gene_ann_path)
    protein_gene = gene_ann[gene_ann.gene_type ==
                            'protein_coding'].gene_name.tolist()
    return(protein_gene)


def _advanced(args, config):
    """functional analysis of a given gene based on correlated coding gene, GSEA analysis"""
    info('running on advanced function')

    config['gene'] = args.gene_name
    config['adj_gene'] = args.adjust
    config['refer_normal'] = args.refer_normal
    config["protein_gene"] = _gene_ann(args.gene_annotation)

    if config['gene']:
        exp_mat = _read_in(config)
        config['exp_mat'] = exp_mat
        config["protein_gene"] = list(
            set(config["protein_gene"]) & set(exp_mat.columns.tolist()))
    else:
        info('gene not given')
        sys.exit(1)
    advanced._run_advanced(config)
