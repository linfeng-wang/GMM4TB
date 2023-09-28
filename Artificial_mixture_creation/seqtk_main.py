#%%
from re import sub
from icecream import ic
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
import re
import os
import json
import subprocess
import itertools
from multiprocessing import pool

#%%
#load the corresponding lineage sample files names and mixture ratio
lin1 = ['ERR752118', 'ERR752118',  'ERR752118', 'ERR2864235', 'ERR2864235', 'ERR2503423', 'ERR751749', 'ERR751749', 'ERR751749', 'ERR2864244', 'ERR2864244', 'ERR2503423',
'ERR2864288', 'ERR2864232', 'ERR2864235', 'ERR2864232', 'ERR2864288', 'ERR2864288', 'ERR2864247', 'ERR2503423', 'ERR2864244', 'ERR2503423', 'ERR2864247', 'ERR2864247',
]

lin2 = ['ERR2864288', 'ERR2864232', 'ERR2864235', 'ERR2864232', 'ERR2864288', 'ERR2864288', 'ERR2864247', 'ERR2503423', 'ERR2864244', 'ERR2503423', 'ERR2864247', 'ERR2864247',
'ERR752118', 'ERR752118',  'ERR752118', 'ERR2864235', 'ERR2864235', 'ERR2503423', 'ERR751749', 'ERR751749', 'ERR751749', 'ERR2864244', 'ERR2864244', 'ERR2503423',

]

ratio1 = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
ratio2 = [100, 95, 90, 85, 80, 75, 70, 65, 60, 55, 50]

# %%
#running the mixture generation script script iteratively
subprocess.run('cd /mnt/storage7/lwang/trial_tb_philippines/pipelines/Genomic_data_analysis/strain_analysis/Artificial_mixture_creation', shell=True)
for l1, l2 in zip(lin1, lin2):
    for r1, r2 in zip(ratio1, ratio2):
        # subprocess.call("./seqtk_new.sh {r1} {r2} {l1} {l2}".format(r1=r1, r2=r2, l1=l1, l2=l2, shell=True))
        subprocess.call("./seqtk_new.sh %s %s %s %s" % (r1, r2, l1, l2), shell=True)
        
#%%
