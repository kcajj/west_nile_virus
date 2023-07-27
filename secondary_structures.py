from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

from mutation_rates import get_mutation_rates, get_genetic_code

def characterise_mutation_site(codon, offset, genetic_code, alphabet):
    synonym=0
    original_nuc=codon[offset]
    original_a=genetic_code[codon]
    for letter in alphabet:
        if letter!=original_nuc:
            if offset==0:
                mutated_codon=letter+codon[1:]
            if offset==1:
                mutated_codon=codon[0]+letter+codon[2]
            if offset==2:
                mutated_codon=codon[:2]+letter
            mutated_a=genetic_code[mutated_codon]
            if original_a==mutated_a:
                synonym+=1
    return synonym

#for each base of the cds compute:
#synonym mutation sites
#non synonym mutation sites
#put them into two lists representing the position of the cds
#for each window on the cds compute the ratio of nonsynonymus mutations over the nonsynonymus sites and 
#the ratio of synonymus mutations over synonyms sites, then compute the ratio of the two ratios

def plot_secondary_structures(k,reference,cds_start,cds_end,alphabet):

    res = get_mutation_rates("wnv/results/ancestral_node_data.json",reference,cds_start,cds_end,alphabet)
    synonym_distribution = res["syn_distribution"]
    genetic_code=get_genetic_code(alphabet)

    cds=reference[cds_start: cds_end]
    synonym_sites_distribution=[0 for i in range(len(synonym_distribution))]
    for pos in range(len(synonym_distribution)):
        codon_offset=pos%3
        codon=cds[pos-codon_offset:pos+(3-codon_offset)]
        
        synonym_site=characterise_mutation_site(codon, codon_offset, genetic_code, alphabet)
        synonym_sites_distribution[pos]=synonym_site

    scores=np.zeros([len(synonym_distribution)-k])
    synonyms=np.zeros([len(synonym_distribution)-k])

    for pos in range(len(synonym_distribution)-k):
        score=sum(synonym_distribution[pos:pos+k])
        synonym_sum=sum(synonym_sites_distribution[pos:pos+k])
        #divide by the number of possible mutations in the window
        scores[pos]=score
        synonyms[pos]=synonym_sum
    
    mean=np.mean(scores)
    print('mean of secondary mutation occurrence:', mean)
    sd=np.std(scores)
    print('sd:',sd)

    sites_without_mutations=[]
    for pos in range(len(synonym_distribution)-k):
        score=sum(synonym_distribution[pos:pos+k])
        z=(score-mean)/sd #normalisation of the observation
        if z<-3:#check if the observation is very small
            sites_without_mutations.append((int(pos+(k/2)),score))

    p=[i for i in range(len(synonym_distribution)-k)]
    secondary_structures=plt.figure(figsize=(25,4))
    plt.plot(p, scores)
    plt.plot(p, synonyms)
    plt.legend(['number of synonymous mutations in the tree','number of available synonymous mutations'])
    for low_score in sites_without_mutations:
        plt.plot(low_score[0],low_score[1],'yo')
    
    return {'sec_struct':sites_without_mutations}