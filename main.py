"""
#!/usr/bin/env python
#author          :Rongbin Zheng
#date            :2020-08-29
#python_version  :3.6 
"""
import sub_functions as sub_cmd
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


class FriendlyArgumentParser(argparse.ArgumentParser):
    """
    Override argparse to show full length help information
    """

    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(1)


"""
subfunctions:
1. distribution: the distribution of the given genes in both tumor and normal
2. boxplot: boxplot between tumor vs normal
3. survial: survival analysis using top 25% vs. bottom 24% samples
4. scatter: scatter plot to show the gene expression correlation between given genes
"""


def parse_args(args=None):
    """
    parse input arguments
    """
    parser = FriendlyArgumentParser(description=__doc__)
    parser.add_argument('-V', '--version', action='version',
                        version="%(prog)s (code version 1.0, db version 1.0)")

    sub_parsers = parser.add_subparsers(
        help="sub-command help", dest="sub_command")

    parser_basic = sub_parsers.add_parser("basic", help="run basic functions",
                                          description="basic functions include distribution, boxplot, survival")
    parser_corr = sub_parsers.add_parser("corr", help="run correlation analysis for a given list of genes",
                                         description="correlation analysis and draw result plot using scatter and heatmap")
    parser_fun = sub_parsers.add_parser("advanced", help="explore the function of the given genes",
                                        description="get gene correlation against a gene, and do GSEA base on the correlation rank")
    # common argument
    for p in [parser_basic, parser_corr, parser_fun]:
        p.add_argument("-i", "--input", dest="inputfile", required=False,
                       help="the path of inputfile where rows are samples, columns are genes, and additional columns include project_id (cancer type), sample_type, survival OS, surival time")

        p.add_argument("-dir", "--inputdir", dest="inputdir", required=False,
                       help="this parameter is lower than -i, the path of input folder where saves each cancer exp matrix file which rows are samples, columns are genes, and additional columns include project_id (cancer type), sample_type, survival OS, surival time")

        p.add_argument("-c", "--cancer", dest="cancer_type", required=False, default='all', type=str,
                       help="specify a cancer type, a list of cancer type, or all cancer type, using comma without space to link cancer name")

        p.add_argument("-p", "--prefix", dest="prefix", required=False,
                       help="the prefix of output files")

        p.add_argument("-ga", "--gene_ann", dest="gene_annotation", required=False, default='./static/gene_annotation_GENCODE.v22.csv', type=str,
                       help="gene annotation where 'gene_name' and 'gene_type' is in the column names, csv file")

    # parser_basic arguments
    parser_basic.add_argument("-d", "--distribution", dest="distribution", required=False, default=False, action="store_true",
                              help="draw expression distribution of given genes in given cancer type")

    parser_basic.add_argument("-b", "--boxplot", dest="boxplot", required=False, default=False, action="store_true",
                              help="draw boxplot to show the expression difference between tumor and normal")

    parser_basic.add_argument("-su", "--survival", dest="survival", required=False, default=False, action="store_true",
                              help="do survival analysis using top 25 percent vs. bottom 25 percent samples based on the given genes")

    parser_basic.add_argument("-g", "--gene", dest="gene_name", required=True, type=str,
                              help="draw expression distribution of given genes in given cancer type")

    parser_basic.add_argument("-m", "--metastatic", dest="metastatic", required=False, default=False, action="store_true",
                              help="plot metastatic and primary comprison in boxplot")

    parser_basic.add_argument("-cli", "--clinic", dest="clinic", required=False,
                              help="clinic table that contains stage information for each sample that match to -i or -dir files, if given the pipeline will plot boxplot in different stage of samples")

    # parser_corr arguments
    parser_corr.add_argument("-g1", "--gene1", dest="gene1", required=True, type=str,
                             help="the first gene")

    parser_corr.add_argument("-g2", "--gene2", dest="gene2", required=True, type=str,
                             help="the second gene")

    parser_corr.add_argument("-j", "--adj", dest="adjust", required=False, type=str,
                             help="the partial correlation will be calculated adjusted by a gene, or mean expression of CD8A,CD8B,GZMB,GZMA,PRF1 if -j CTL is set")

    # evaluate the function of given genes using GSEA
    parser_fun.add_argument('-g', '--gene', dest='gene_name', required=True,
                            help="the gene name that you want to explore the function")

    parser_fun.add_argument("-j", "--adj", dest="adjust", required=False,
                            help="the partial correlation will be calculated adjusted by a gene, or CTL level (mean expression of CD8A,CD8B,GZMB,GZMA,PRF1)")

    parser_fun.add_argument("-rn", "--refer_normal", dest="refer_normal", required=False, action='store_true',
                            help="gained positive / negative correlation relative to normal sample if normal is available")

    return (parser.parse_args(args), parser)


def info(mesg):
    os.system('echo "++++ %s"' % mesg)


def main(args=None):
    args, parser = parse_args(args)
    if not args.sub_command:
        os.system('python main.py -h')
        sys.exit(1)
    config = {}  # add argument in dict format
    if not args.inputfile and not args.inputdir:
        info('Please set -i or -d')
        sys.exit(1)
    config['sub_cmd'] = args.sub_command
    config['fpath'] = args.inputfile
    config['dpath'] = args.inputdir
    # all cancer type if not specified
    config['cancer'] = args.cancer_type if args.cancer_type else 'all'
    # give NA as prefix if not specified
    config['prefix'] = args.prefix if args.prefix else 'NA'
    # =============== basic module ===============
    if args.sub_command == 'basic':
        sub_cmd._basic(args, config)

    # =============== correlation analysis module ===============
    if args.sub_command == 'corr':
        sub_cmd._corr(args, config)

    # =============== advanced functional analysis ===============
    if args.sub_command == 'advanced':
        sub_cmd._advanced(args, config)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr("User interrupt:) \n")
