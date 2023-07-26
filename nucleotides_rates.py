import json
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def create_mut_matrix(alphabet):
    m={}
    for i in alphabet:
        for j in alphabet:
            if not i==j:
                m[i+j]=0
    return m


def get_nuc_usage(reference,cds_start,cds_end):
    nuc_counts_s={'A':0,'C':0,'T':0,'G':0}
    #counts the nucleotides occurrence i
    for pos, nuc in enumerate(reference):
        if pos>=cds_start and pos<cds_end: #we do not consider the upper limit because our sequence starts from zero, the last nucleotide of the cds is 10394
            if (pos-cds_start)%3==2:
                nuc_counts_s[nuc]+=1
    nuc_usage_s={'A':0,'C':0,'T':0,'G':0}
    total_positions_s=0
    for count in nuc_counts_s.values():
        total_positions_s+=count
    for nuc in nuc_usage_s:
        nuc_usage_s[nuc]=nuc_counts_s[nuc]/total_positions_s
    
    nuc_counts_n={'A':0,'C':0,'T':0,'G':0}
    #counts the nucleotides occurrence i
    for pos, nuc in enumerate(reference):
        if pos>=cds_start and pos<cds_end: #we do not consider the upper limit because our sequence starts from zero, the last nucleotide of the cds is 10394
            if not((pos-cds_start)%3==2):
                nuc_counts_n[nuc]+=1
    nuc_usage_n={'A':0,'C':0,'T':0,'G':0}
    total_positions_n=0
    for count in nuc_counts_n.values():
        total_positions_n+=count
    for nuc in nuc_usage_n:
        nuc_usage_n[nuc]=nuc_counts_n[nuc]/total_positions_n

    return {'nuc_usage_syn':nuc_usage_s,'nuc_usage_nonsyn':nuc_usage_n}
