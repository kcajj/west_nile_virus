import json
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

from nucleotides_rates import get_nuc_usage
from mutation_rates import get_genetic_code, get_mutation_rates
from eigen_decomposition import plot_GTR
from secondary_structures import plot_secondary_structures
from KaKs_ratio import plot_dnds_ratio

alphabet='ACTG'
cds_start=96 #the sequence annotation is in base 1, here the sequence starts from 0, so the first codon is 96-97-98
cds_end=10395
reference=SeqIO.read('wnv/config/reference.gb','gb').seq
k=50
tree_node_data="wnv/results/ancestral_node_data.json"

nuc_empirical_sequences=get_nuc_usage(reference,cds_start,cds_end,alphabet)

plot_GTR(reference,cds_start,cds_end,alphabet,tree_node_data)

secondary=plot_secondary_structures(k,reference,cds_start,cds_end,alphabet,tree_node_data)
print(secondary['sec_struct'])
print('the secondary structure seems to extend from', 3475+96, 'to', 3523+96)

dnds=plot_dnds_ratio(k,reference,cds_start,cds_end,alphabet,tree_node_data)
print(dnds['avg_ratio'])

plt.show()