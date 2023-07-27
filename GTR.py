import numpy as np
import matplotlib.pyplot as plt

from nucleotides_rates import create_nuc_count_matrix, get_nuc_usage
from mutation_rates import create_mut_matrix, get_mutation_rates

def plot_GTR(reference,cds_start,cds_end,alphabet,tree_node_data):

    nuc_usage=get_nuc_usage(reference,cds_start,cds_end,alphabet)
    res = get_mutation_rates(tree_node_data,reference,cds_start,cds_end,alphabet)
    synonym_counts = res["syn_counts"]
    total_synonym_mut = res['tot_syn_mut']

    dictQ=create_mut_matrix(alphabet)
    for mut in dictQ:
        dictQ[mut]=(synonym_counts[mut]/total_synonym_mut)
    
    diagonal_rates=create_nuc_count_matrix(alphabet)
    for mut in dictQ:
        diagonal_rates[mut[0]]+=-dictQ[mut]

    #adding diagonal frequences
    for letter in diagonal_rates:
        dictQ[letter+letter]=diagonal_rates[letter]

    #from dictionary to matrix
    print('matrix Q:')
    Q=np.zeros([4,4])
    for i in range(len(alphabet)):
        for j in range(len(alphabet)):
            Q[i,j]=(dictQ[alphabet[i]+alphabet[j]])
    Q=Q.T
    print(Q)

    control_eq_p=[[],[]]
    for i in nuc_usage['nuc_usage_nonsyn']:
        control_eq_p[0].append(nuc_usage['nuc_usage_syn'][i])
        control_eq_p[1].append(nuc_usage['nuc_usage_nonsyn'][i])

    print('initial conditions')
    initial=np.array([1,0,0,0])
    print(initial)

    print('eigen values')
    e_val=np.linalg.eigvals(Q)
    print(e_val)

    #######
    #EIGEN DECOMPOSITION
    #######

    revect=np.linalg.eig(Q)[1]
    levect=np.linalg.eig(Q.T)[1]
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
    GTR=plt.figure()
    plt.plot(t, [probabilities_in_function_of_time(ti,a,e_val,revect,) for ti in t])
    plt.legend(['A','C','T','G'])