from pgmpy.models import BayesianNetwork
import argparse
import os
import argparse
import pandas as pd
import sys
import networkx as nx
import matplotlib.pyplot as plt
import random
import seaborn as sns
from typing import List, Tuple
import numpy as np
from prob_agg import compute_union_graph, graph_agony
import itertools
from pgmpy.models import BayesianNetwork
from pgmpy.estimators import HillClimbSearch,BDeuScore, K2Score, BicScore
import ioutils as io
import tronco_utils as tu
import networkutils as nu
from pgmpy.models import BayesianNetwork
import pandas as pd
import itertools
import os
# authors @ abhinav tamaskar and james bannon
import cvxopt as cx
import networkx as nx
import numpy as np
import random
import math
import sys


class ExtendedBN(BayesianNetwork):
    """
        Extends BayesianNetwork Class to include a numbered networkx graph
        and a label dictionary
    """
    def __init__(self,edge_list:str,mutations:str):
        pass

    def make_pgm(edge_list:str,mutations:str)->BayesianNetwork:
        mutation_data = pd.read_csv(mutations,index_col=0).reset_index(drop=True)
        nodes = list(mutation_data.columns)
        model = BayesianNetwork()
        model.add_nodes_from(nodes)
        f = open(edge_list, "r")
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
        for edge in lines[3:]:
            edge=edge.split()
            model.add_edge(edge[0],edge[1])
        model.fit(mutation_data)
        return model


    def get_event_probability(model,evidence:dict)-> float:
        log_sum = 0
        # print("event probs")
        nodes = model.nodes()
        for node in nodes:
            cpd = model.get_cpds(node)
            vals = cpd.get_values()
            parents = model.get_parents(node)
            if len(parents)>=1:
                cpd_=cpd.reduce([(node,evidence[node]) for node in parents],inplace=False)
                vals = cpd_.get_values()
            log_sum += np.log(vals[evidence[node]][0])
        return np.exp(log_sum)


    def compute_tv_vector(model:BayesianNetwork)->np.array:
        # print(model.nodes)
        if len(model.nodes) > 20:
            all_events = sample_binary_events(model.nodes)
        else:
            all_events = pd.DataFrame(data=np.array(list(map(list, itertools.product([0, 1],
            repeat=len(model.nodes))))),columns = model.nodes)

        results = []
        for i,r in all_events.iterrows():
            res = get_event_probability(model,r.to_dict())
            results.append(res)
        return np.array(results)

    def get_event_probability(model,evidence:dict)-> float:
        log_sum = 0
        # print("event probs")
        nodes = model.nodes()
        for node in nodes:
            cpd = model.get_cpds(node)
            vals = cpd.get_values()
            parents = model.get_parents(node)
            if len(parents)>=1:
                cpd_=cpd.reduce([(node,evidence[node]) for node in parents],inplace=False)
                vals = cpd_.get_values()
            log_sum += np.log(vals[evidence[node]][0])
        return np.exp(log_sum)

    #TODO: make sure it allows passsing in `all' as an argument
    def fit_bayesnet(mutations:pd.DataFrame,fname:str) ->None:
        print("***fitting bn***")
        print(fname)

        model = HillClimbSearch(mutations)
        model = model.estimate(scoring_method=BicScore(mutations))
        model = BayesianNetwork(model)
        model.fit(mutations)
        nx.draw(model, with_labels=True)
        plt.savefig(fname)
        plt.close()
        return model

    def read_edge_list_to_nx_graph(edge_list:str)->nx.DiGraph:
        D=nx.DiGraph()
        f = open(edge_list, "r")
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
        num_nodes = lines[2]
        for j in range(0,int(num_nodes)):
            D.add_node(j)
        num_edges = lines[1]
        edges = lines[3:]
        for edge in edges:
            # print(edge)
            edge= edge.split()
            D.add_edge(int(edge[0]),int(edge[1]))
        return D
