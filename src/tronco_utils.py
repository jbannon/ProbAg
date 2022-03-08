from typing import Tuple, List
import pandas as pd
import os
import argparse
import sys
import networkx as nx
from pgmpy.models import BayesianNetwork
from pgmpy.estimators import HillClimbSearch,BDeuScore, K2Score, BicScore
import ioutils as io
from BayesianNetwork import ExtendedBN

def consolidate_cancer(mutation_info:pd.DataFrame)->Tuple[List[str],List[str]]:

    fname = "../swapspace/temp_cons.csv"
    mutation_info.to_csv(fname)
    os.system("Rscript consolidate.r {fn}".format(fn=fname))
    os.remove(fname)

    if os.path.exists("../swapspace/keep.txt"):
        keep = io.get_list("../swapspace/keep.txt")
        toss = io.get_list("../swapspace/remove.txt")
        os.remove("../swapspace/keep.txt")
        os.remove("../swapspace/remove.txt")
    else:
        keep = []
        toss = []

    return keep,toss


def fit_sbcn(data:pd.DataFrame,
    exp_name:str,
    group_name:str,
    cancer:str,
    group:str,
    filter:str,
    cutoff:int) -> BayesianNetwork:

    f = open("../swapspace/sbcn_info.txt","w")
    f.writelines([x+"\n" for x in [exp_name,group_name,cancer,filter,str(cutoff)]])
    f.close()

    fname = "../swapspace/mutations.csv"
    data.to_csv(fname)
    os.system("Rscript fit_sbcn.r {fn}".format(fn=fname))
    SBCN = bu.make_pgm("../swapspace/edge_list.txt",fname)
    NX = bu.read_edge_list_to_nx_graph("../swapspace/edge_list_numeric.txt")
    os.remove(fname)
    os.remove("../swapspace/sbcn_info.txt")
    os.remove("../swapspace/edge_list.txt")
    os.remove("../swapspace/edge_list_numeric.txt")
    return SBCN,NX
    #TODO: move edge list
    #os.system("mv ../swapspace/edge_list.txt ../results/{e}/{g}/")
