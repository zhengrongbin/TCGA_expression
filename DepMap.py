import os
import sys
import math
import re
import pandas as pd
import numpy as np
import collections
import pickle as pk
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

import pingouin as pg
# import statsmodels.api as sm

# from gprofiler import GProfiler
plt.rcParams.update(plt.rcParamsDefault)

# plt.style.use('seaborn-ticks')
sns.set(style='ticks')
rc = {'axes.labelpad': 15, "axes.labelsize": 16,
      "figure.titleweight": "bold"}
# # mpl.style.use('ggplot')
# rc={"axes.labelsize": 16, "xtick.labelsize": 12, "ytick.labelsize": 12,
#     "figure.titleweight":"bold", #"font.size":14,
#     "figure.figsize":(5.5,4.5), "font.weight":"regular", "legend.fontsize":10,
#    'axes.labelpad':10,
# }
plt.rcParams.update(**rc)

gene_name = sys.argv[1] #official gene symbol in DepMap data matrix

# crispr screening result
depmap_screening = pd.read_csv(
    './static/DepMap_Screen_result.csv.gz', compression='gzip', index_col=0)
# sample annotation
sample_info_path = './static/sample_info.csv'
sample_info = pd.read_csv(sample_info_path, index_col=0)
# make model matrix for correcting by tissue type
tmp = sample_info[['sample_collection_site']]
tmp = tmp[tmp.index.isin(depmap_screening.index)]
tmp['stat'] = 1
tmp = tmp.pivot(columns='sample_collection_site', values='stat')
tmp = tmp.fillna(0)
tmp = tmp.iloc[:, 1:].astype('int')

# ============= plot distribution
fig, ax = plt.subplots()
ax.hist(depmap_screening[gene_name], density=True)
sns.despine(ax=ax, offset=5, trim=True)
ax.set(xlabel='Score', ylabel='Frequency')
ax.set_title(label='DepMap Screening in CCLE cell lines',
             fontdict={'fontweight': 'bold', 'fontsize': 14})
# plt.show()
plt.tight_layout()
fig.savefig('_distribution.pdf')
plt.close()

# calculate partial correlation of all genes besides given gene, adjusted by tissue types
corr_res = {}
for g in depmap_screening.columns.tolist():
    if g == gene_name:
        continue
    df = pd.concat([depmap_screening[[gene_name, g]].astype(
        'float'), tmp.astype('int')], axis=1).dropna()
    try:
        corr_res[g] = pg.partial_corr(df.astype('float'),
                                      x=g, y=gene_name, covar=tmp.columns.tolist()).loc['pearson', ['r', 'p-val']]
    except:
        continue

corr_res = pd.DataFrame.from_dict(corr_res, orient='index')
corr_res.columns = ['pcorr', 'ppval']
corr_res = corr_res.sort_values('pcorr', ascending=False)
corr_res.to_csv('DepMap_screen_partial_corr_adjByTissueType.csv')
