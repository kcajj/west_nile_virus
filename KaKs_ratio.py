import matplotlib.pyplot as plt

from mutation_rates import get_mutation_rates, get_genetic_code

def characterise_mutation_site(codon, offset, genetic_code, alphabet):
    synonym=0
    non_synonym=0
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
            else:
                non_synonym+=1
    return (synonym/(synonym+non_synonym), non_synonym/(synonym+non_synonym))

#for each base of the cds compute:
#synonym mutation sites
#non synonym mutation sites
#put them into two lists representing the position of the cds
#for each window on the cds compute the ratio of nonsynonymus mutations over the nonsynonymus sites and 
#the ratio of synonymus mutations over synonyms sites, then compute the ratio of the two ratios

def plot_dnds_ratio(k,reference,cds_start,cds_end,alphabet):

    res=get_mutation_rates("wnv/results/ancestral_node_data.json",reference,cds_start,cds_end,alphabet)
    synonym_distribution=res['syn_distribution']
    non_synonym_distribution=res['nonsyn_distribution']
    genetic_code=get_genetic_code(alphabet)

    cds=reference[cds_start: cds_end]
    synonym_sites_distribution=[0 for i in synonym_distribution]
    non_synonym_sites_distribution=[0 for i in synonym_distribution]
    for pos in range(len(synonym_distribution)):
        codon_offset=pos%3
        codon=cds[pos-codon_offset:pos+(3-codon_offset)]
        
        synonym_site, non_synonym_site=characterise_mutation_site(codon, codon_offset, genetic_code, alphabet)
        synonym_sites_distribution[pos]+=synonym_site
        non_synonym_sites_distribution[pos]+=non_synonym_site

    ratios=[]
    for pos in range(0,len(synonym_distribution)-k,3):
        try:
            ds=sum(synonym_distribution[pos:pos+k])/sum(synonym_sites_distribution[pos:pos+k])
            dn=sum(non_synonym_distribution[pos:pos+k])/sum(non_synonym_sites_distribution[pos:pos+k])
            ratio=dn/ds
        except ZeroDivisionError:
            ratio=0
        ratios.append(ratio)

    average_ratio=sum(ratios)/len(ratios)

    p=[i for i in range(0,len(synonym_distribution)-k,3)]
    kd_ks_ratio=plt.figure(figsize=(25,4))
    plt.plot(p, ratios)
    plt.legend(['dn/ds ratio'])

    return {'syn_sites_distr':synonym_sites_distribution,'nonsyn_sites_distr':non_synonym_sites_distribution,'avg_ratio':average_ratio}
