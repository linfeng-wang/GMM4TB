#must input mixed infection cells
#%%
import argparse
from array import array
from statistics import mode
import numpy as np #numpy version 1.22 needed
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
# import seaborn as sns
# import matplotlib.pyplot as plt
import pandas as pd
# from icecream import ic
from uuid import uuid4
import sys
from pathlib import Path

import tb_profiler 
import gmm_model_multi as gmm_model
import pprint


#%% testing
# json_file = '../strain_analysis/test_data/ERR6634978-ERR6635032-3070.results.json' #file used for targeting and error checking
# vcf_file = '/mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis/test_data/ERR6634978-ERR6635032-3070.vcf.gz' #file used creating the model

# # vcf_file = '/mnt/storage7/jody/tb_ena/per_sample/ERR221637.gatk.vcf.gz'
# # json_file='/mnt/storage7/jody/tb_ena/tbprofiler/latest/results/ERR221637.results.json'

# vcf_file = '/mnt/storage7/lwang/Projects/Philipine_tb_report/vcf/freebayesVCF_new/ERR6635386.vcf.gz'
# json_file = '/mnt/storage7/jody/philippines/projects/all_philippines/tbprofiler_new/results/ERR6635386.results.json'

# graph_option = False
# output_path = '/mnt/storage7/lwang/Projects/Philipine_tb_report/test'
# output_name = 'test3'

#%% #Input commands don't run in test
output_path = None
multi = False
vcf_file = None
fq1_file = None
fq2_file = None
json_file = None
vcf_input = False

###
parser = argparse.ArgumentParser(description='Mixed_infection_GMM',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-vcf", "--vcf", help='VCF (gatk) file')
parser.add_argument("-v", "--use_vcf_input", help='Use VCF as input file to generate TB-profiler json file', action='store_true')
parser.add_argument("-json", "--json", help='TB-profiler output json file')
parser.add_argument("-fq1", "--fastq1", help='Fastq forward read file')
parser.add_argument("-fq2", "--fastq2", help='Fastq backward read file')
# parser.add_argument("-json", "--json", help='tb-profiler output json file')
parser.add_argument("-m", "--multi_infection_detection", help='consider 2+ strain mixed infection as well', action='store_true')
parser.add_argument("-g", "--graph", help='alternative snp frequency histogram', action='store_true',default=False)
parser.add_argument("-o", "--output_name", help='output name')
parser.add_argument("-op", "--output_path", help='output path')

args = parser.parse_args()

vcf_input = args.use_vcf_input
graph_option = args.graph
json_file1 = args.json
vcf_file = args.vcf
output_name = args.output_name
output_path = args.output_path
multi = args.multi_infection_detection
fq1_file = args.fastq1
fq2_file = args.fastq2
#%%
subprocess.run(f"mkdir -p {output_path}/temp", shell=True)

if vcf_file != None and vcf_input == True and json_file1 == None:
    # subprocess.run(f"mkdir -p {output_path}/temp", shell=True)
    subprocess.run(f"tb-profiler profile --vcf {vcf_file} --dir {output_path}/temp/ -p {output_name}", shell=True)
    json_file = f'{output_path}/temp/results/{output_name}.results.json'

    # subprocess.run(f"tb-profiler profile --vcf {vcf_file} -p {output_name} --dir {output_path}/temp/", shell=True)
elif json_file1 != None:
    json_file = json_file1
# if fq1_file != None and fq2_file != None:
#     subprocess.run("python variant_calling_from_fastq.py -1 {fq1_file} -2 {fq2_file} -r $REFGENOME -o temp/{output_name}", shell=True)
#     subprocess.run(f"tb-profiler profile --vcf temp/{output_name}.vcf.gz -p temp/{output_name}",  shell=True)
#%%    
# print(json_file)
#%%
tb_pred_result, output_status = tb_profiler.tb_pred(json_file)

#%%
if output_status == 1 and multi == True:
    print("=======================")
    print(f"Programme continued, 2+ strain mixture found in {vcf_file}")
    print("=======================")
    
elif output_status == 1 and multi == False:
    print("=======================")
    sys.exit(f"Programme stoped, contamination! in {vcf_file} (2+ mixture found)")
    print("=======================")
    
elif output_status == 2:
    sys.exit(f"Programme stoped, no mixed infection in {vcf_file}")
    
elif vcf_file == None and (fq1_file == None or fq2_file == None) and json_file == None:
    sys.exit(f"Programme stoped, no input file indicated")

else:
    print("***********************")
    print(f"Programme continued, mixed infection detected in {vcf_file}")
    print("***********************")

dr_dict = tb_profiler.tb_dr(json_file)
# subprocess.run(f"rm -r {output_path}/temp/", shell=True)

#%%
gmm_pred_result, model = gmm_model.model_pred(vcf_file, 
                                            tail_cutoff = list(tb_pred_result.values())[1], 
                                            graph = False, 
                                            output_path=output_path, 
                                            mix_num=len(tb_pred_result))

print('mix_num:',len(tb_pred_result))
print('gmm_pred_result:',gmm_pred_result)

if sum(gmm_pred_result) < 0.9: #adding threshold if the sum of the fraction lower than 0.9, then rejected
    sys.exit(f"Programme stoped, low prediction confidence in {vcf_file}")

gmm_pred_result, model = gmm_model.model_pred(vcf_file, 
                                            tail_cutoff = list(tb_pred_result.values())[1], 
                                            graph = graph_option, 
                                            output_path=output_path, 
                                            mix_num=len(tb_pred_result))

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

#%%
strain_dict = {}

for x in range(len(gmm_pred_result)):
    strain_dict[f"s{x}"] = {}
unknown = []

for x in list(strain_dict.keys()):
    strain_dict[x]["DR_pred"] = []
# assigning resistance to strain first this strain arrangement is the same as strain the output from gmm
for element in dr_dict: #key is freqs value is dr
    prob = model.predict_proba(np.array(element["freq"]).reshape(-1,1))
    prob = prob[0] #the output from predict_proba is a list of a list
    ind_max_prob = np.argmax(prob)
    dict_ = {"prediction_confidence" : prob[ind_max_prob], 
                "info" : element}
    # the below if statement is to check if the prediction confidence is higher than 0.99, if so then it will be assigned to the strain with the highest prediction confidence
    if np.array(element["freq"]) > 0.99:
        for k,v in strain_dict.items():
            strain_dict[k]["DR_pred"].append(dict_)
    else: 
        # strain_dict[list(strain_dict.keys())[ind_max_prob]]["DR_pred"] = []
        strain_dict[list(strain_dict.keys())[ind_max_prob]]["DR_pred"].append(dict_)
        
        if all(elem == prob[0] for elem in prob):
            dict_ = {"prediction_confidence" : prob, 
                    "info" : element}
            unknown.append(dict_)

#%% effort trying to sort the the dr pred in order of prediction_confidence - doesn't work now
# for k, v in strain_dict.items():
#     v['DR_pred'] = sorted(v["DR_pred"], key=lambda d: d[0]['prediction_confidence'], reverse=True) 
    #print(v)
    # print(v['DR_pred'])
    # x["DR_pred"] = sorted(x["DR_pred"], key=lambda d: d['prediction_confidence']) 

m_mean = np.concatenate(model.means_)
model_covariances_ = np.concatenate(np.concatenate(model.covariances_))
m_covar = np.sqrt(model_covariances_/2) #standard deviation

for x in range(len(gmm_pred_result)):
    strain_dict[list(strain_dict.keys())[x]]["gmm_pred"] = m_mean[x]
    strain_dict[list(strain_dict.keys())[x]]["gmm_UB"] = m_mean[x]+m_covar[x]
    strain_dict[list(strain_dict.keys())[x]]["gmm_LB"] = m_mean[x]-m_covar[x]

#%% failed attempt to order by gmm pred
# strain_dict = dict(sorted(strain_dict.items(), key=lambda item: item["gmm_pred"], reverse=True))
# strains[0] = dict(sorted(strains[0].items(), key=lambda item: item[1], reverse=True))

#So the order of out pred_prob output should be the same as model_pred results, hence which adjust the ordering of these two together
#tb prediction results has already by adjusted to have the bigger value infront so it doesn't need to be adjusted this way.

#%% using bubble sort to order by gmm pred
n = len(strain_dict)

for i in range(0,n): #bubble sort to get the lineage list from high to low gmm predicted proportion wise
    for j in range(i+1, n):
        a = (strain_dict[f"s{i}"]["gmm_pred"])
        b = (strain_dict[f"s{j}"]["gmm_pred"])
        if a < b:
            strain_dict[f"s{i}"], strain_dict[f"s{j}"] = strain_dict[f"s{j}"], strain_dict[f"s{i}"]

for x in range(len(gmm_pred_result)):
    strain_dict[list(strain_dict.keys())[x]]["lin"] = list(tb_pred_result.keys())[x]
    strain_dict[list(strain_dict.keys())[x]]["tb_pred"] = list(tb_pred_result.values())[x]

#%% 
# #adjusting the order of the output presented
strain_dict1 = {}
for x in range(len(gmm_pred_result)):
    strain_dict1[f"s{x}"] = {}
    strain_dict1[f"s{x}"]["lin"] = strain_dict[f"s{x}"]["lin"]
    strain_dict1[f"s{x}"]["tb_pred"] = strain_dict[f"s{x}"]["tb_pred"]
    strain_dict1[f"s{x}"]["gmm_pred"] = strain_dict[f"s{x}"]["gmm_pred"]
    strain_dict1[f"s{x}"]["gmm_UB"] = strain_dict[f"s{x}"]["gmm_UB"]
    strain_dict1[f"s{x}"]["gmm_LB"] = strain_dict[f"s{x}"]["gmm_LB"]
    strain_dict1[f"s{x}"]["DR_pred"] = strain_dict[f"s{x}"]["DR_pred"]

#%%
dr_output = {}
dr_output["sample_name"] = output_name
dr_output["lineage"] = []
dr_output["tb_pred"] = []
dr_output["gmm_pred"] = []
for x in strain_dict1.values():
    dr_output["lineage"].append(x['lin'])
    dr_output["tb_pred"].append(x['tb_pred'])
    dr_output["gmm_pred"].append(x['gmm_pred'])

dr_output["mixture"] = strain_dict1

dr_output["Unknown50-50DR"] = unknown
dr_output["gmm_tb-profiler_MSE"]= mse


desired_order = ['sample_name', 'lineage', 'tb_pred', 'gmm_pred', 'gmm_tb-profiler_MSE', 'mixture', 'Unknown50-50DR', ]

# Create a new dictionary with the desired order
dr_output = {key: dr_output[key] for key in desired_order}
#Outputting result into below format
    # dr_output = {"lineage": [
    #                 {
    #                 "lin" : list(tb_pred_result.keys())[0],
    #                 "tb_pred" : list(tb_pred_result.values())[0],
    #                 "gmm_pred" :  m_mean[0],
    #                 "gmm_pred_UB" : strain0_bound[0],
    #                 "gmm_pred_LB" : strain0_bound[1],
    #                 "resistance_pred" : strains_0
    #                 },
    #                 {
    #                 "lin" : list(tb_pred_result.keys())[1],
    #                 "tb_pred" : list(tb_pred_result.values())[1],
    #                 "gmm_pred" :  m_mean[1],
    #                 "gmm_pred_UB" : strain1_bound[0],
    #                 "gmm_pred_LB" : strain1_bound[1],
    #                 "resistance_pred" : strains_1
    #                 }
    # ],
    #             "Unknown50-50DR" : unknown,
    #             "gmm_tb-profiler_MSE": mse
    # }
    
#%%
# name = vcf_file.split('/')[-1].split('.vcf')[0]
name = output_name

if output_path == None:
    cwd = os.path.dirname(__file__)
    Path(f"{cwd}/mixture_results/profiles").mkdir(parents=True, exist_ok=True)
    output_file = "".join([f"{cwd}/mixture_results/profiles","/",name, ".mix.json"])

else:
    Path(f"{output_path}/mixture_results/profiles").mkdir(parents=True, exist_ok=True)
    output_file = "".join([f"{output_path}/mixture_results/profiles","/",name, ".mix.json"])

with open(output_file, 'w') as f:
    json.dump(dr_output, f, indent=2)

#%%
print("=======================================================================")
print(name)
print("*************************")
pprint.pprint(dr_output)
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

#vcf file
# python /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable/gmm_model_tbp/main_multi.py -vcf /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis/test_data/ERR6634978-ERR6635032-3070.vcf.gz -g -op /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/Executable_eval -o vcf_input_test


#fasq