#must input mixed infection cells
#%%
import argparse
from statistics import mode
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
from scipy.stats import kurtosis, skew
import re
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from icecream import ic
from uuid import uuid4
import sys
from pathlib import Path

#%%
import tb_profiler 
import gmm_model

#%% testing
# json_file = '../strain_analysis/test_data/ERR6634978-ERR6635032-3070.results.json' #file used for targeting and error checking
# vcf_file = '../strain_analysis/test_data/ERR6634978-ERR6635032-3070.vcf.gz' #file used creating the model
# graph_option = False
# output_path = None


#%%
output_path = None

parser = argparse.ArgumentParser(description='Mixed_infection_GMM',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-vcf", "--vcf", help='VCF (gatk) file')
parser.add_argument("-json", "--json", help='tb-profiler output json file')
parser.add_argument("-g", "--graph", help='alternative snp frequency histogram', action='store_true')
parser.add_argument("-o", "--output", help='output path')



args = parser.parse_args()

graph_option = args.graph
json_file = args.json
vcf_file = args.vcf
output_path = args.output

#%%
tb_pred_result, failed = tb_profiler.tb_pred(json_file)

if failed == 1:
    sys.exit(f"Programme stoped, contamination! in {vcf_file}")
elif failed == 2:
    sys.exit(f"Programme stoped, no mixed infection in {vcf_file}")
else:
    print("***********************")
    print(f"Programme continued, mixed infection detected in {vcf_file}")
    print("***********************")

dr_dict = tb_profiler.tb_dr(json_file)

gmm_pred_result, model = gmm_model.model_pred(vcf_file, tail_cutoff = list(tb_pred_result.values())[1], graph = graph_option, output_path=output_path)

if sum(gmm_pred_result) < 0.9: #adding threshold if the sum of the fraction lower than 0.9, then rejected
    sys.exit(f"Programme stoped, low prediction confidence in {vcf_file}")

mse = gmm_model.mse_cal(tb_pred_result, gmm_pred_result)

#%%###########################################################################
# strains = list(tb_pred_result.keys())
# strains[0] = {}
# strains[1] = {}
# unknown = {}

# for key, value in dr_dict.items(): #key is freqs value is dr
#     prob = model.predict_proba(np.array(float(key)).reshape(-1,1))
#     prob = prob[0] #the output from predict_proba is a list of a list
#     if prob[0] > prob[1]:
#         strains[0][value] = prob[0]
#     elif prob[0] < prob[1]:
#         strains[1][value] = prob[1]
#     else:
#         unknown[value] = prob

# strains = list(tb_pred_result.keys())
# strains[0] = []
# strains[1] = []
strains_0 = []
strains_1 = []
unknown = []

for element in dr_dict: #key is freqs value is dr
    prob = model.predict_proba(np.array(element["freq"]).reshape(-1,1))
    prob = prob[0] #the output from predict_proba is a list of a list
    if np.array(element["freq"]) > 0.99: # if frequency is very high then the same mutation may be in both lineages
        dict_ = {"prediction_confidence" : np.max(prob), 
                "info" : element}
        strains_0.append(dict_)
        
        dict_ = {"prediction_confidence" : np.max(prob), 
                "info" : element}
        strains_1.append(dict_)
        
    elif prob[0] > prob[1]:
        dict_ = {"prediction_confidence" : prob[0], 
                "info" : element}
        strains_0.append(dict_)

    elif prob[0] < prob[1]:
        dict_ = {"prediction_confidence" : prob[1], 
                "info" : element}
        strains_1.append(dict_)
    
    else:
        dict_ = {"prediction_confidence" : prob, 
                "info" : element}
        unknown.append(dict_)

strains_0 = sorted(strains_0, key=lambda d: d['prediction_confidence']) 
strains_1 = sorted(strains_1, key=lambda d: d['prediction_confidence']) 

# strains[0] = dict(sorted(strains[0].items(), key=lambda item: item[1], reverse=True))
# strains[1] = dict(sorted(strains[1].items(), key=lambda item: item[1], reverse=True))


m_mean = np.concatenate(model.means_)
model_covariances_ = np.concatenate(np.concatenate(model.covariances_))
m_covar = np.sqrt(model_covariances_/2 ) #standard deviation

strain0_bound = [m_mean[0]+m_covar[0], m_mean[0]-m_covar[0]] #1 std interval upper and lower bound for model predicted proportion for both strains
strain1_bound = [m_mean[1]+m_covar[1], m_mean[1]-m_covar[1]]

#So the order of out pred_prob output should be the same as model_pred results, hence which adjust the ordering of these two together
#tb prediction results has already by adjusted to have the bigger value infront so it doesn't need to be adjusted this way.
#%%
if m_mean[0] >= m_mean[1]:

    dr_output = {"lineage": [
                    {
                    "lin" : list(tb_pred_result.keys())[0],
                    "tb_pred" : list(tb_pred_result.values())[0],
                    "gmm_pred" :  m_mean[0],
                    "gmm_pred_UB" : strain0_bound[0],
                    "gmm_pred_LB" : strain0_bound[1],
                    "resistance_pred" : strains_0
                    },
                    {
                    "lin" : list(tb_pred_result.keys())[1],
                    "tb_pred" : list(tb_pred_result.values())[1],
                    "gmm_pred" :  m_mean[1],
                    "gmm_pred_UB" : strain1_bound[0],
                    "gmm_pred_LB" : strain1_bound[1],
                    "resistance_pred" : strains_1
                    }
    ],
                "Unknown50-50DR" : unknown,
                "gmm_tb-profiler_MSE": mse
    }

else:

    dr_output = {"lineage": [
                    {
                    "lin" : list(tb_pred_result.keys())[0],
                    "tb_pred" : list(tb_pred_result.values())[0],
                    "gmm_pred" :  m_mean[1],
                    "gmm_pred_UB" : strain1_bound[0],
                    "gmm_pred_LB" : strain1_bound[1],
                    "resistance_pred" : strains_1
                    },
                    {
                    "lin" : list(tb_pred_result.keys())[1],
                    "tb_pred" : list(tb_pred_result.values())[1],
                    "gmm_pred" :  m_mean[0],
                    "gmm_pred_UB" : strain0_bound[0],
                    "gmm_pred_LB" : strain0_bound[1],
                    "resistance_pred" : strains_0
                    }
    ],
                "Unknown50-50DR" : unknown,
                "gmm_tb-profiler_MSE": mse
    }


#%%
# print(dr_output)        
#%%
            # list(tb_pred_result.keys())[1]:strains[1]}

name = vcf_file.split('/')[-1].split('.vcf')[0]

if output_path == None:
    cwd = os.path.dirname(__file__)
    Path(f"{cwd}/mixture_results/profiles").mkdir(parents=True, exist_ok=True)
    output_file = "".join([f"{cwd}/mixture_results/profiles","/",name, ".mix.json"])

else:
    Path(f"{output_path}/mixture_results/profiles").mkdir(parents=True, exist_ok=True)
    output_file = "".join([f"{output_path}/mixture_results/profiles","/",name, ".mix.json"])

with open(output_file, 'w') as f:
    json.dump(dr_output, f, indent=2)

print("=======================================================================")
print(name)
print(dr_output)
print("=======================================================================")


# def tb_profiler_dr(json_file):
#     for var in json_results['dr_variants']:
#         cluster = assign_variant_to_distrib(gm, var['freq'], cutoff=0)
#         if cluster[0] == 0:
#             var['probs']= cluster
#             strain0.append(var)
#         elif cluster[0] == 1:
#             var['probs']= cluster
#             strain1.append(var)
#         else:
#             var['probs']= cluster[1:]
#             strainU.append(var)

#     strain_0 = [(v['drugs'][0]['drug']) for v in strain0]
#     strain_1 = [(v['drugs'][0]['drug']) for v in strain1]

# # %%
# def assign_variant_to_distrib(gm,freq,cutoff=0.95):
#     probs = gm.predict_proba([[freq]])
#     pred = gm.predict([[freq]])
#     if probs[0][pred][0]>cutoff: #put the prediction probability through cutoff threshold
#         return pred[0], probs[0][pred][0]
#     else:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
#         return None, pred[0], probs[0][pred][0]

# %%
#command to run
# python /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/main.py -vcf /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis/test_data/ERR6634978-ERR6635032-3070.vcf.gz -json /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis/test_data/ERR6634978-ERR6635032-3070.results.json -g -o /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable_eval


# python /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/main.py -vcf /mnt/storage7//jody/tb_ena/per_sample/ERR2864229.gatk.vcf.gz -json /mnt/storage7/jody/tb_ena/tbprofiler/latest/results/ERR2864229.results.json -g -o /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable_eval
# python /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/main.py -vcf /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Closeup/vcf/ERR5897084.gatk.vcf.gz -json /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Closeup/results/ERR5897084.results.json -g -o /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable_eval

# python /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/main.py -vcf /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Closeup/vcf/ERR4829977.gatk.vcf.gz -json /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Closeup/results/ERR4829977.results.json -g -o /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable_eval
