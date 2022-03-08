from collections import namedtuple
from typing import NamedTuple, Union, List
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys


class MetaData(NamedTuple):
    """ container for cancer server info"""
    tissue:str
    subset:str
    is_breast:bool

class CancerDataServer:
    """
        Class to perform cancer data fetching, preprocessing, and plotting utils
        keeps a set of binarized data frames vaf threshold.
        subsetting/filting genes must be done by building an external list and then
        calling subset mutations
    """
    def __init__(self,
    tissue:str,
    subset:str = "all",
    ) -> None:

        self.metadata = MetaData(*[tissue,subset,tissue.lower()=='brca'])
        self.mutations = pd.read_csv("../data/{c}/preprocessed/mutation_vaf.csv".format(c=tissue))
        self.get_patient_ids()
        self.vaf_to_binary_df = {}


    def fit_from_vaf_thresh_list(self,thresholds:List[float]):
        for t in tqdm(thresholds,desc="thresholds"):
            self.binarize(t)

    def binarize(self,thresh:float):
        if thresh in self.vaf_to_binary_df.keys():
            return self.vaf_to_binary_df[thresh]
        else:
            temp = self.mutations
            temp['Above Thresh'] = np.where(self.mutations['VAF']>thresh,1,0)
            events = list(pd.unique(temp['Hugo_Symbol']))
            event_to_idx  = {}
            c=0
            for e in events:
                event_to_idx[e]=c
                c+=1
            rows = {}
            for i in tqdm(pd.unique(temp['Tumor_Sample_Barcode']),desc='Patients'):
                patient_df = temp[(temp['Tumor_Sample_Barcode']==i) & (temp['Above Thresh']==1)]
                gene_names = list(pd.unique(patient_df['Hugo_Symbol']))
                row = np.zeros(len(events),dtype=int)
                for m in gene_names:
                    row[event_to_idx[m]]=1
                rows[i]=row
            res = pd.DataFrame.from_dict(rows,orient='index',columns=events)
            self.clean(res)
            self.vaf_to_binary_df[thresh] = res

    def get_patient_ids(self):
        if self.metadata.subset.lower()=="all":
            self.patient_ids = self.mutations['Tumor_Sample_Barcode'].to_list()
        else:
            pass
            """TODO: implement logic to get the right subset patient ids"""

    def subset_mutations(self,vaf_key:float,genes:List[str])->pd.DataFrame:
        if vaf_key in self.vaf_to_binary_df.keys():
            df = self.vaf_to_binary_df[key]
            return df[set(df.columns).intersection(set(genes))]

    def get_mutations(self)->pd.DataFrame:
        return self.mutations

    def get_binarized_mutations(self,vaf_key:float)->pd.DataFrame:
        if vaf_key in self.vaf_to_binary_df:
            return self.vaf_to_binary_df[vaf_key]


    def clean(self,data:pd.DataFrame)->pd.DataFrame:
        """ assumes binarized df """
        colsums = data.sum()
        events_in_all_samples = colsums[colsums==data.shape[0]].index
        events_that_occur = colsums.iloc[colsums.to_numpy().nonzero()].index

        data = data[events_that_occur]

        if len(events_in_all_samples)>0:
            data.drop(columns=events_in_all_samples,inplace=True)
        return data




if __name__ == '__main__':
    g = CancerDataServer('brca','all')
    g.fit_from_vaf_thresh_list(VAF_THRESHOLDS)
    # for t in VAF_THRESHOLDS:
    #     tt = g.get_binarized_mutations(t)
    #     print(tt.head())
