import itertools
import numpy as np
from tqdm import tqdm
from typing import List, Union, Dict, Tuple
from CancerDataServer import CancerDataServer
import tronco_utils as tu
import matplotlib.pyplot as plt
import pandas as pd
import os
import argparse
import sys
import networkx as nx
from pgmpy.models import BayesianNetwork
from pgmpy.estimators import HillClimbSearch,BDeuScore, K2Score, BicScore
import ioutils as io

VAF_MIN = 0
VAF_MAX = 0.5
VAF_STEP = 0.1
VAF_THRESHOLDS = list(np.arange(VAF_MIN,VAF_MAX+VAF_STEP,VAF_STEP))
CANCERS = ['brca','ov','gbm']

FILTER_DICT = {'t1':'t1_genes.txt','t2':'t2_genes.txt',
    'hm':'hallmark_genes.txt','t1h':'t1_hallmarks.txt','t2h':'t2_hallmarks.txt'}


def process_cancer_pair(
    c1_mutations:pd.DataFrame,
    c2_mutations:pd.DataFrame,
    cutoff:int
    )->Tuple[pd.DataFrame,pd.DataFrame]:
    """
    Take in a pair of data frames read from specific cancer mutation files and
    process them. Specifically for each cancer we:

    1) clean (remove absent and ever-present genes)
    2) keep only genes/events which have measured in both and by (1) observed in both
    3) consolidate each of them separately using TRONCO

    """

    common_genes = list(set(c1_mutations.columns).intersection(set(c2_mutations.columns)))
    c1_mutations = c1_mutations[common_genes]
    c2_mutations = c2_mutations[common_genes]
    master_df = c1_mutations.append(c2_mutations,ignore_index=True)

    colsums = master_df.sum()
    master_df = master_df[colsums.iloc[np.lexsort([colsums.index, -colsums.values])].index[:cutoff]]

    c1_mutations = c1_mutations[list(master_df.columns)]
    c2_mutations = c2_mutations[list(master_df.columns)]

    keep = []
    toss = []
    k1,t1 = tu.consolidate_cancer(c1_mutations)
    keep.extend(k1)
    toss.extend(t1)

    k2,t2 = tu.consolidate_cancer(c2_mutations)
    keep.extend(k2)
    toss.extend(t2)

    final_kept_columns = list(master_df.columns)
    for x in toss:
        final_kept_columns.remove(x)
    master_df = master_df[final_kept_columns]

    return c1_mutations[list(master_df.columns)],c2_mutations[list(master_df.columns)]



def get_filters(file_table:Dict[str,str])->Dict[str,Union[List[str],str]]:
    filter_lists = {}
    for key in file_table.keys():
        filter_lists[key] = io.get_list("../data/filters/{k}".format(k=file_table[key]))
    filter_list['none']=None
    return filter_list




def build_data_servers(cancers:List[str],vaf_thresholds:List[float]) -> Dict[str,CancerDataServer]:
    print("*************************************************\n*\t                            \t\t*\n*\tBuilding Cancer Data Servers\t\t*\n*\n*\n*************************************************")
    CDS_table = {}

    for cancer in cancers:
        print(cancer)
        CDS_table[cancer]= CancerDataServer(cancer,"all")
        CDS_table[cancer].fit_from_vaf_thresh_list(VAF_THRESHOLDS)
    return CDS_table

def create_experiment_directories(
    exp_name:str,
    group_name:str,
    cancer1:str,
    group1:str,
    cancer2:str,
    group2:str,
    filter:str,
    cutoff:int
    ) -> None:

    exp_name = "{c1}_v_{c2}".format(c1=cancer1,c2=cancer2)
    group_name = "{g1}_v_{g2}".format(g1=group1,g2=group2)
    os.makedirs("../figs/{en}/{gn}/{f}/{k}/{c}/".format(en=exp_name,gn=group_name,f=str(filter),k=str(cutoff),c=cancer1),exist_ok=True)
    os.makedirs("../figs/{en}/{gn}/{f}/{k}/{c}/".format(en=exp_name,gn=group_name,f=str(filter),k=str(cutoff),c=cancer2),exist_ok=True)

    os.makedirs("../results/{en}/{gn}/{f}/{k}/{c}/".format(en=exp_name,gn=group_name,f=str(filter),k=str(cutoff),c=cancer1),exist_ok=True)
    os.makedirs("../results/{en}/{gn}/{f}/{k}/{c}/".format(en=exp_name,gn=group_name,f=str(filter),k=str(cutoff),c=cancer2),exist_ok=True)





def make_figures(
    cancers: List[str],
    vaf_thresholds:List[float],
    number_genes_to_keep: int = 15,
    filter:str= 'none'
    )->None:


    # ToDo: Cache once experiment is done

    cancer_pairs = itertools.combinations(cancers,2)
    CDS_table = build_data_servers(cancers,vaf_thresholds)
    for cancer1, cancer2 in itertools.combinations(cancers,2):
        CS1 = CDS_table[cancer1]
        CS2 = CDS_table[cancer2]
        for thresh in VAF_THRESHOLDS:
            df1 = CS1.get_binarized_mutations(thresh)
            df2 = CS2.get_binarized_mutations(thresh)
            df1, df2 = process_cancer_pair(df1,df2,number_genes_to_keep)
            print(df1.head())






if __name__ == '__main__':
    make_figures(CANCERS,VAF_THRESHOLDS)
