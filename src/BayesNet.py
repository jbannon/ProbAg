import os
import pandas as pd
import sys
import networkx as nx
from typing import List, Tuple, Union
import numpy as np
import itertools
from pgmpy.models import BayesianNetwork
from pgmpy.estimators import HillClimbSearch,BDeuScore, K2Score, BicScore
from CancerDataServer import CancerDataServer

class ExtendedBN:
    """
        Extends BayesianNetwork Class to include a
        numbered networkx graph and a label dictionary
    """
    def __init__(self):
        self.fit=False

    def get_tv_vector(self)->np.array:
        if self.TV_vect is None:
            self.TV_vect = self.compute_tv_vector()
        return self.TV_vect

    def fit_as_sbcn(self,data:pd.DataFrame)->None:
        fname = "../swapspace/mutations.csv"
        data.to_csv(fname)

        os.system("Rscript fit_sbcn.r {fn}".format(fn=fname))

        self.make_name_maps("../swapspace/node_numbering.txt")

        nodes = list(data.columns)
        model = BayesianNetwork()
        model.add_nodes_from(nodes)
        model.fit(data)
        self.model = model
        self.read_edge_list_to_nx_graph("../swapspace/edge_list_numeric.txt")

        os.remove(fname)
        os.remove("../swapspace/edge_list.txt")
        os.remove("../swapspace/edge_list_numeric.txt")


    def make_name_maps(self,fname:str):
        gene_2_idx, idx_2_gene = {},{}

        with open("../swapspace/node_numbering.txt","r") as f:
            lines = f.readlines()
            lines = [line.rstrip() for line in lines]
            for line in lines:
                gene, idx = line.split()
                gene_2_idx[gene] = idx
                idx_2_gene[idx] = gene
        self.gene_2_idx = gene_2_idx
        self.idx_2_gene = idx_2_gene


    @classmethod
    def fit_as_bn(cls,mutation_data:pd.DataFrame)->None:
        model = HillClimbSearch(mutation_data)
        model = model.estimate(scoring_method=BicScore(mutation_data))
        model = BayesianNetwork(model)
        model.fit(mutation_data)
        self.model = model
        self.fit = True

    def _echo_maps(self):
        print(self.gene_2_idx)
        print(self.idx_2_gene)

    def compute_tv_vector(self)->np.array:
        if len(self.model.nodes) > 20:
            all_events = sample_binary_events(model.nodes)
        else:
            all_events = pd.DataFrame(data=np.array(list(map(list, itertools.product([0, 1],
            repeat=len(model.nodes))))),columns = model.nodes)

        results = []
        for i,r in all_events.iterrows():
            res = self.get_event_probability(model,r.to_dict())
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

    def get_nx_graph(self)->nx.DiGraph:
        return self.nx_graph

    def read_edge_list_to_nx_graph(self,edge_list:str)->None:
        D=nx.DiGraph()
        print(edge_list)
        with open(edge_list, "r") as f:
            lines = f.readlines()
            lines = [line.rstrip() for line in lines]
        num_nodes = lines[2]
        for j in range(0,int(num_nodes)):
            D.add_node(j)
        num_edges = lines[1]
        edges = lines[3:]
        for edge in edges:
            edge= edge.split()
            D.add_edge(int(edge[0]),int(edge[1]))
        self.nx_graph = D

if __name__ == '__main__':
    cancer = 'brca'
    vts = [0.1,0.2]
    S = CancerDataServer(cancer)
    S.fit_from_vaf_thresh_list(vts)
    for vt in vts:
        df = S.get_binarized_mutations(vt)
        df = df.iloc[:,:10]
        BN = ExtendedBN()
        BN.fit_as_sbcn(df)
        G = BN.get_nx_graph()
        print(G)
        BN._echo_maps()
