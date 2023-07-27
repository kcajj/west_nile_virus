import json
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

from nucleotides_rates import get_nuc_usage
from mutation_rates import create_mut_matrix, get_mutation_rates

alphabet='ACTG'
cds_start=96 #the sequence annotation is in base 1, here the sequence starts from 0, so the first codon is 96-97-98
cds_end=10395
reference=SeqIO.read('wnv/config/reference.gb','gb').seq

nuc_usage=get_nuc_usage(reference,cds_start,cds_end)
print(nuc_usage)
res = get_mutation_rates("wnv/results/ancestral_node_data.json",reference,cds_start,cds_end,alphabet)
synonym_counts = res["syn_counts"]
total_synonym_mut = res['tot_syn_mut']

Q=create_mut_matrix(alphabet)
for mut in Q:
    Q[mut]=(synonym_counts[mut]/total_synonym_mut)
print('just mutation frequences')
print(Q)
diagonal_rates={'A':0,'C':0,'T':0,'G':0}
for mut in Q:
    diagonal_rates[mut[0]]+=-Q[mut]

print('also diagonal frequences')
for letter in diagonal_rates:
    Q[letter+letter]=diagonal_rates[letter]
print(Q)

print('matrix Q:')
matrix=np.empty([4,4])
for i in range(len(alphabet)):
    row=[]
    for j in range(len(alphabet)):
        matrix[i,j]=(Q[alphabet[i]+alphabet[j]])
matrix=matrix.T
print(matrix)

print('this should be equal to the equilibrium')
control_eq_p=[[],[]]
for i in nuc_usage['nuc_usage_nonsyn']:
    control_eq_p[0].append(nuc_usage['nuc_usage_syn'][i])
    control_eq_p[1].append(nuc_usage['nuc_usage_nonsyn'][i])
print(control_eq_p)

print('initial conditions')
initial=np.array([1,0,0,0])
print(initial)

print('eigen values')
e_val=np.linalg.eigvals(matrix)
print(e_val)

#######
#EIGEN DECOMPOSITION
#######
revect=np.linalg.eig(matrix)[1]
levect=np.linalg.eig(matrix.T)[1]
print('rev')
print(revect[:,0]/sum(revect[:,0]))

def eigen_projection(p,lv,rv):
    return [np.dot(w,p)/np.dot(w,v) for w,v in zip(lv.T,rv.T)]

print('eigen decomposition of Q')
a=eigen_projection(initial,levect,revect)
print(a)

#rec=np.sum([a[i]*v for i,v in enumerate(revect.T)], axis=0)
#print(np.round(rec,3),initial)

def probabilities_in_function_of_time(t,a,ev,rv):
    return np.sum([a[i]*np.exp(ev[i]*t)*v for i,v in enumerate (rv.T)],axis=0)

print('equilibrium probabilities')
print(probabilities_in_function_of_time(100, a, e_val, revect))
print('empirical probabilities (first synonymous then non synonymous positions)')
print(control_eq_p[0])
print(control_eq_p[1])

t=np.linspace(0,100,100)
plt.plot(t, [probabilities_in_function_of_time(ti,a,e_val,revect,) for ti in t])
plt.legend(['A','C','T','G'])
plt.show()

