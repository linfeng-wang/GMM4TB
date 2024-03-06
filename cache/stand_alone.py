#%%
import argparse
from ntpath import join
import uuid
import numpy as np
# import pathogenprofiler as pp
from sklearn.mixture import GaussianMixture
from sklearn.metrics import mean_squared_error, confusion_matrix, f1_score
from sklearn.metrics import make_scorer
from sklearn.model_selection import cross_val_score

# import fastq2matrix as fm
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
from scipy.stats import norm
import subprocess
# from scipy.stats import kurtodsdsis, skew
import re
import os
# import seaborn as sns
# import matplotlib.pyplot as plt
import pandas as pd
# from cb91visuals import *
# import pipe
# from icecream import ic
from pathlib import Path

#%%
#this model outputs the strain info depending if one value is still larger than the threshold after substracting the other prob value from it - not as good as not using it

# #%% test
# vcf_file = '/mnt/storage7/jody/tb_ena/per_sample/ERR221662.gatk.vcf.gz' #file used creating the model
# vcf_file = '/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk/freebayesVCF/ERR6634978-ERR6635032-3070.vcf.gz'
# vcf_file = '/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk/freebayesVCF/ERR6634978-ERR6635032-5050.vcf.gz'
# vcf_file = '/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk/freebayesVCF/ERR6634978-ERR6635032-1000.vcf.gz'
# vcf_file = '/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk/freebayesVCF/ERR6634978-ERR6635032-9505.vcf.gz'

#%%
parser = argparse.ArgumentParser(description='Stand_alone_mix_predictor',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input_csv", help='input csv file with two columns, one is the position or SNP name, the other is the alternative allele frequency')
parser.add_argument("-c", "--mix_count", help='User suggested count of mixes (default=2)', default=2, type=int)
parser.add_argument("-o", "--output_path", help='file path to save the output')
parser.add_argument("-g", "--graphics", help='alternative snp frequency histogram', action='store_true', default=False)

args = parser.parse_args()

input = args.input_csv
mix_count = args.mix_count
output_path = args.output_path
graph = args.graphics
#%%
#test inputs
# input = '/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/gmm_model_tbp/cache/stand_alone_test.txt'
# output_path = 'stand_alone_test_out.csv'
# mix_count = 2

#%%
df = pd.read_csv(input, header=0)
a = np.array(df.iloc[:,1]).reshape(-1, 1)
gm = GaussianMixture(n_components=mix_count, random_state=0).fit(a)
cluster_no = []

for x in a:
    prob = gm.predict_proba([x])
    cluster_no.append(np.argmax(prob)+1)

print('='*20)
print('User suggested mix count:', mix_count)
print('Lineage fractions:', np.concatenate(gm.means_).tolist())
print('Output file:', output_path)
print('='*20)

df['cluster_no'] = cluster_no
df.to_csv(output_path, index=False)
    # if graph:
    
model_covariances =[]
model_means_ = np.concatenate(gm.means_)
model_covariances_ = np.concatenate(np.concatenate(gm.covariances_))
model_covariances_ = np.sqrt(model_covariances_/2 )

for i, x in enumerate(model_covariances_):
    model_covariances.append([model_means_[i]+model_covariances_[i], model_means_[i]-model_covariances_[i]])

# plt.savefig(output_file, format="pdf", bbox_inches="tight")
#%%
if graph:
    flat_freqs = list(np.concatenate(a))
    fig = go.Figure()
    fig.add_trace(go.Histogram(x = flat_freqs,nbinsx=round(len(flat_freqs)/10)))

    fig.update_layout(
        title="Alternative allele GMM distribution",
        xaxis_title="Alternative SNP freqency<br><sup>Highlights show interval +/- 1STD around the mean as determined by gmm</sup>",
        yaxis_title="Count",
        font=dict(
            family="Courier New, monospace",
            size=18,
            color="RebeccaPurple"
        )
    )

    for x in model_covariances:
        fig.add_vrect(x0=x[0], x1=x[1], row="all", col=1,
            annotation_text="", annotation_position="top left",
            fillcolor="green", opacity=0.25, line_width=0)
        
    name = output_path.split('/')[-1].split('.csv')[0]
    output_file = "".join([name, "_gmm_fig.png"])
    fig.write_image(output_file)

#test run
# python stand_alone.py -i /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/gmm_model_tbp/cache/stand_alone_test.txt -c 2 -o stand_alone_test_out.csv -g
