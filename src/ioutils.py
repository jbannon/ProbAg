from typing import List, Tuple
import argparse
import pandas as pd


def get_list(fname:str)->List[str]:
    f = open(fname,"r")
    lines = f.readlines()
    lines = [line.rstrip() for line in lines]
    f.close()
    return lines

def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c1','--cancer1',type=str)
    parser.add_argument('-g1','--group1',type=str)
    parser.add_argument('-c2','--cancer2',type=str)
    parser.add_argument('-g2','--group2',type=str)
    parser.add_argument('-f','--filter',type=str)
    parser.add_argument('-k','--cutoff',type=int)
    return parser



def make_dirs(cancer1:str,cancer2:str,filter:str,cutoff:int):
    exp_name = "{c1}_v_{c2}".format(c1=cancer1,c2=cancer2)
    os.makedirs("./figs/{en}/{c1}/{g1}/{f}/{k}".format(en=exp_name,c1=cancer1,g1=group1,f=filter,k=str(cutoff)),exist_ok=True)
    os.makedirs("./figs/{en}/{c2}/{g2}/{f}/{k}".format(en=exp_name,c2=cancer2,g2=group2,f=filter,k=str(cutoff)),exist_ok=True)
    #os.makedirs("./figs/{c1}_{c2}/{c1}/{g}/{f}/{k}".format(c1=cancer1,c2=cancer2,g=group2,f=filter,k=str(cutoff)),exist_ok=True)
    os.makedirs("./results/{en}/{f}/{k}".format(en=exp_name,f=filter,k=str(cutoff)),exist_ok=True)
