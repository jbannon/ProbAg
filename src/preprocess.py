"""
Code to pre-store VAF computation
"""
import pandas as pd
import numpy as np
import os
import sys

def compute_and_write_vaf(
    cancer:str
    )->None:
    dataset = pd.read_table("../data/{c}/data_mutations.txt".format(c=cancer))
    dataset = dataset[['Hugo_Symbol',
        'Variant_Classification',
        'Variant_Type',
        't_ref_count',
        't_alt_count',
        'Tumor_Sample_Barcode',
        'Mutation_Status']]
    print(cancer)
    print(len(pd.unique(dataset['Tumor_Sample_Barcode'])))
    dataset['VAF'] = dataset['t_alt_count']*1.0/(dataset['t_ref_count']+dataset['t_alt_count'])
    dataset.to_csv("../data/{c}/preprocessed/mutation_vaf.csv".format(c=cancer),index=False)

def BRCA_subtyping(dataset:pd.DataFrame):
    dataset=dataset.dropna(subset=['SUBTYPE'])
    count_table = dataset['SUBTYPE'].value_counts()
    count_table.to_csv("../data/brca/preprocessed/subtype_counts.csv",index=False)

    for subtype in pd.unique(dataset['SUBTYPE']):
            temp = dataset[dataset['SUBTYPE']==subtype]
            st_label = subtype[5:]
            temp.to_csv("../data/brca/preprocessed/clinical_{st}.csv".format(st=st_label),index=False)


def process_clinical_ids(
    cancer:str,
    nquantiles:int = 4
    ) -> None:

    os.makedirs("../data/{c}/preprocessed/known_pfs/".format(c=cancer),exist_ok=True)
    os.makedirs("../data/{c}/preprocessed/all_pfs/".format(c=cancer),exist_ok=True)

    sample = pd.read_table("../data/{c}/data_clinical_sample.txt".format(c=cancer),skiprows=4)
    patient = pd.read_table("../data/{c}/data_clinical_patient.txt".format(c=cancer),skiprows=4)

    if cancer =='brca':
            BRCA_subtyping(patient)


    patient_sample_map = sample[['PATIENT_ID','SAMPLE_ID']]
    patient_sample_map.to_csv("../data/{c}/preprocessed/id_map.csv".format(c=cancer))


    patient = patient[['PATIENT_ID','SUBTYPE','PFS_STATUS','PFS_MONTHS']]

    known_pfs = patient[patient["PFS_STATUS"]!='0:CENSORED']
    labels = list(np.arange(nquantiles)+1)
    labels.reverse()
    known_pfs['PFS_QTILE']=pd.qcut(known_pfs['PFS_MONTHS'],q=nquantiles,labels=labels)
    for label in labels:
        temp = known_pfs[known_pfs['PFS_QTILE']==label]
        temp.to_csv("../data/{cancer}/preprocessed/known_pfs/clinical_q{q}.csv".format(cancer=cancer,q=label),index=False)
    count_table =known_pfs['PFS_QTILE'].value_counts()
    count_table.to_csv("../data/{cancer}/preprocessed/known_pfs/quantile_counts.csv".format(cancer=cancer),index=False)


    patient['PFS_QTILE']=pd.qcut(patient['PFS_MONTHS'],q=nquantiles,labels=labels)
    for label in labels:
        temp = patient[patient['PFS_QTILE']==label]
        temp.to_csv("../data/{cancer}/preprocessed/all_pfs/clinical_q{q}.csv".format(cancer=cancer,q=label),index=False)
    count_table =patient['PFS_QTILE'].value_counts()
    count_table.to_csv("../data/{cancer}/preprocessed/all_pfs/quantile_counts.csv".format(cancer=cancer),index=False)



if __name__ == '__main__':

    for cancer in ['brca','gbm','kirc','lgg','luad','lusc','ov','prad','thca','ucec']:
        print(cancer)
        os.makedirs("../data/{c}/preprocessed/".format(c=cancer),exist_ok=True)
        compute_and_write_vaf(cancer)
        process_clinical_ids(cancer)
