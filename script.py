import json
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

#creates an empty matrix 4X4
def create_mut_matrix(alphabet):
    m={}
    for i in alphabet:
        for j in alphabet:
            if not i==j:
                m[i+j]=0
    return m

alphabet='ACTG'
cds_start=96
cds_end=10395
synonym_counts=create_mut_matrix(alphabet)

children_of_root=[]
#finds the children of the root and adds them to a list
with open('wnv/auspice/wnv.json') as auspice:
    tree_data=json.load(auspice)
    tree=tree_data['tree']
    children=tree['children']
    for child in children:
        children_of_root.append(child['name'])

#counts the mutations, the root mutations are counted half in one direciton and half in the other
with open('wnv/results/ancestral_node_data.json') as augur:
    node_data = json.load(augur)
    nodes = node_data['nodes']
    for node in nodes:
        mutations=nodes[node]['muts']
        for mutation in mutations:

            from_nt=mutation[0]
            to_nt=mutation[-1]
            position=int(mutation[1:-1])

            if position>=cds_start and position<cds_end:
                if (position-cds_start)%3==0:
                    if node in children_of_root:
                        synonym_counts[from_nt+to_nt]+=1/2
                        synonym_counts[to_nt+from_nt]+=1/2
                    else:
                        synonym_counts[from_nt+to_nt]+=1

print(synonym_counts)

reference=SeqIO.read('wnv/data/reference_sequence.fna','fasta')
nuc_counts={'A':0,'C':0,'T':0,'G':0}
#counts the nucleotides occurrence i
for pos, nuc in enumerate(reference.seq):
    if pos>=cds_start and pos<cds_end:
        if (pos-cds_start)%3==2:
            nuc_counts[nuc]+=1

nuc_usage={'A':0,'C':0,'T':0,'G':0}
total_positions=0
for i in nuc_counts.values():
    total_positions+=i
for nuc in nuc_usage:
    nuc_usage[nuc]=nuc_counts[nuc]/total_positions

print(nuc_usage)

ratios_matrix=create_mut_matrix(alphabet)
total_synonym_mut=0
for i in synonym_counts.values():
    total_synonym_mut+=i
for mut in ratios_matrix:
    ratios_matrix[mut]=(synonym_counts[mut]/total_synonym_mut)*nuc_usage[mut[0]]

print(ratios_matrix)

print()
print()
print()


#my try of computing the equilibrium probabilities

Q=create_mut_matrix(alphabet)
for mut in Q:
    Q[mut]=(synonym_counts[mut]/total_synonym_mut)
print('just mutation frequences')
print(Q)
diagonal_probabilities={'A':0,'C':0,'T':0,'G':0}
for mut in Q:
    diagonal_probabilities[mut[0]]+=-Q[mut]

print('also diagonal frequences')
for letter in diagonal_probabilities:
    Q[letter+letter]=diagonal_probabilities[letter]
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
control_eq_p=[]
for i in nuc_usage:
    control_eq_p.append(nuc_usage[i])
print(control_eq_p)

print('initial conditions')
initial=np.array([1,0,0,0])
print(initial)

print('eigen values')
e_val=np.linalg.eigvals(matrix)
print(e_val)

'''
print('eigen vectors')
e_vect=np.linalg.eig(matrix)[1]
print(e_vect)

t0=0

p0=np.empty(4)
for i in range(len(e_vect[0])):
    p=0
    for j in range(len(e_val)):
        p+=e_vect[j,i]*(np.exp(e_val[j]*t0))*initial[j]
    p0[i]=p

print('p0')
print(p0)

def final_function(t,p0,e_val,e_vect):
    return np.sum([p0[i]*np.exp(e_val[i]*t)*v for i,v in enumerate(e_vect.T)],axis=0)

print('equilibrium probabilities')
print(final_function(100, p0, e_val, e_vect))
#t=np.linspace(0,100,100)
#plt.plot(t, [final_function(ti,p0,e_val,e_vect) for ti in t])
#plt.show()

#i_matrix=np.array([[e_val[0],0,0,0],[0,e_val[1],0,0],[0,0,e_val[2],0],[0,0,0,e_val[3]]])
#for i in range(len(e_val)):
#    i_matrix[i,i]=e_val[i]

#def fun(t,e_vect,i_matrix):
#    poft=np.dot(np.dot(np.linalg.inv(e_vect),np.exp(i_matrix*t)),e_vect)
#    solutions=np.linalg.solve(poft,[0,0,0,0])
#    return solutions
#print(fun(100, e_vect, i_matrix))


'''

#######
#EIGEN DECOMPOSITION
#######
print()
print()
print()

revect=np.linalg.eig(matrix)[1]
levect=np.linalg.eig(matrix.Tx)[1]
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
print('my probabilities')
print(control_eq_p)
t=np.linspace(0,100,100)
plt.plot(t, [probabilities_in_function_of_time(ti,a,e_val,revect,) for ti in t])
plt.legend(['A','C','T','G'])
plt.show()

