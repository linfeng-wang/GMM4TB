#%%
import argparse
from ntpath import join
import uuid
import numpy as np
# import pathogenprofiler as pp
from sklearn.mixture import GaussianMixture
from sklearn.metrics import mean_squared_error, confusion_matrix, f1_score
# import fastq2matrix as fm
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
from scipy.stats import norm
import subprocess
# from scipy.stats import kurtodsdsis, skew
import re
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
# from cb91visuals import *
import gc
import pipe
from icecream import ic
from uuid import uuid4
from pathlib import Path


#%%
#this model outputs the strain info depending if one value is still larger than the threshold after substracting the other prob value from it - not as good as not using it

#%% test
vcf_file = '../strain_analysis/test_data/ERR6634978-ERR6635032-2080.vcf.gz' #file used creating the model
# vcf_file = '/mnt/storage7/jody/tb_ena/per_sample/ERR221662.gatk.vcf.gz' #file used creating the model
vcf_file = '/mnt/storage7/lwang/trial_tb_philippines/data/processed/seqtk/freebayesVCF/ERR6634978-ERR6635032-3070.vcf.gz'
#%%
def model_pred(vcf_file, tail_cutoff=0.07, graph = False, output_path = None, mix_num = 3):
    cwd = os.path.dirname(__file__) #this is used to get the folder location of the script so that new_exculsion file can be accessed

    if not os.path.exists(f"{cwd}/temp"):
        os.mkdir(f"{cwd}/temp")

    uuid_ = uuid4() #create random naming so that the script can be run in parallel
    uuid_file = "".join([str(uuid_), ".csv"])


    with open(f"{cwd}/temp/{uuid_file}", 'w') as f:
        subprocess.run(f"bcftools view -c 1 -m2 -M2 -T ^{cwd}/new_exclusion.bed %s | bcftools query -f '%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD\\n]'" % vcf_file, shell=True, stdout=f, text=True)
    pos = []
    freqs = []
    scatter = []
    with open(f"{cwd}/temp/{uuid_file}", 'r') as f:
        for line in f: #get the relevant info including position, alt/ref snp count info

            row = line.strip().split()
            ads = [int(x) for x in row[4].split(",")]
            afs = [x/(sum(ads)+1e-14) for x in ads]
            if afs[1]>1-tail_cutoff or afs[1]<tail_cutoff: #apply cut off value
                continue
            pos.append(int(row[0]))
            freqs.append([afs[1]])
            scatter.append(ads)

        # freqs = [[0.7],[0.6],[0.4]]    
        gm = GaussianMixture(n_components=mix_num, random_state=0).fit(np.array(freqs).reshape(-1, 1))
        mu0 = gm.means_[1][0]
        mu1 = gm.means_[0][0]
        
        scatter = np.array(scatter)
        
        labels = gm.predict(freqs)
        mu0_ = len(labels[labels==0])/len(labels)
        mu1_ = len(labels[labels==1])/len(labels)
        
        prediction_prob = gm.predict_proba(freqs)
        
        
        anchored_array = np.concatenate((scatter, labels.reshape(-1,1)), axis=1) 
        lin1_ = anchored_array[anchored_array[:,2]==1]
        lin4_ = anchored_array[anchored_array[:,2]==0]

    if graph:

        # flat_freqs = list(np.concatenate(freqs))
        # plt.hist(flat_freqs)
        # plt.title("Alternative allele GMM distribution")
        # plt. xlabel("Alternative SNP freqency")
        # plt. ylabel("Count")
        model_means_ = np.concatenate(gm.means_)
        model_covariances_ = np.concatenate(np.concatenate(gm.covariances_))
        model_covariances_ = np.sqrt(model_covariances_/2 )
        strain1_bound = [model_means_[0]+model_covariances_[0], model_means_[0]-model_covariances_[0]] #1 std interval upper and lower bound for model predicted proportion for both strains
        strain2_bound = [model_means_[1]+model_covariances_[1], model_means_[1]-model_covariances_[1]]

        # plt.savefig(output_file, format="pdf", bbox_inches="tight")

        flat_freqs = list(np.concatenate(freqs))
        fig = go.Figure()
        fig.add_trace(go.Histogram(x = flat_freqs))
    
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
    
        fig.add_vrect(x0=strain1_bound[0], x1=strain1_bound[1], row="all", col=1,
            annotation_text="", annotation_position="top left",
            fillcolor="green", opacity=0.25, line_width=0)

        fig.add_vrect(x0=strain2_bound[0], x1=strain2_bound[1], row="all", col=1,
            annotation_text="", annotation_position="top left",
            fillcolor="green", opacity=0.25, line_width=0)

        name = vcf_file.split('/')[-1].split('.vcf')[0]
        output_file = "".join([name, ".gmm_fig.png"])

        if output_path == None:
            Path(f"{cwd}/mixture_results/images").mkdir(parents=True, exist_ok=True)
            fig.write_image(f"{cwd}/mixture_results/images/{output_file}")
        else:
            Path(f"{output_path}/mixture_results/images").mkdir(parents=True, exist_ok=True)
            fig.write_image(f"{output_path}/mixture_results/images/{output_file}")


        # fig.show()

    # os.remove(os.path.join("/__pycache__", uuid_file))
    os.remove(f"{cwd}/temp/{uuid_file}")

    if mu0 > mu1:
        return [mu0, mu1], gm #make sure the descending order corresponds to the output of tb_profiler.tb_pred()
    else:
        return [mu1, mu0], gm

#Test
# model_pred(vcf_file, tail_cutoff=0, graph = False)


# %%
def mse_cal(tb_prof, gmm_result):
    return mean_squared_error(list(tb_prof.values()), gmm_result)


#%% test
# bcftools view -c 1 -m2 -M2 -T ^/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/new_exclusion.bed ERR2503421.g.vcf.gz | bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT\\t%AD\\n]' | head -1
#bcftools view -c 1 -m2 -M2 -T DRR014508.g.vcf.gz

# bcftools view -c 1 -m2 -M2 -T ^/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/new_exclusion.bed ../strain_analysis/test_data/ERR6634978-ERR6635032-3070.vcf.gz | bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT\\t%AD\\n]'

# bcftools view -c 1 -m2 -M2 -T ^/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/new_exclusion.bed /mnt/storage7//jody/tb_ena/per_sample/ERR2864229.gatk.vcf.gz | bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT\\t%AD\\n]' | head -1
