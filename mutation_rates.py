import json
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup

from nucleotides_rates import create_mut_matrix

alphabet='ACTG'
cds_start=96 #the sequence annotation is in base 1, here the sequence starts from 0, so the first codon is 96-97-98
cds_end=10395
reference=SeqIO.read('wnv/config/reference.gb','gb').seq

def get_genetic_code():
    genetic_code={}
    with open('genetic-code.html') as page:
        gencodepage=BeautifulSoup(page,'html.parser')
    table_rows=gencodepage.find_all('tr')
    for row in table_rows:
        cells=row.find_all('td')
        if cells != []: #remove some empty rows
            codon=cells[0].get_text()
            if codon[0] in alphabet: #ensure that the row is valid
                aminoacid=cells[2].get_text()
                genetic_code[codon]=aminoacid
    return genetic_code

def is_synonym(from_codon, to_codon, genetic_code):
    return genetic_code[from_codon]==genetic_code[to_codon]

def get_mutation_rates(mutation_json_name):
    
    children_of_root=[]
    #finds the children of the root and adds them to a list
    with open('wnv/auspice/wnv.json') as auspice:
        tree_data=json.load(auspice)
        tree=tree_data['tree']
        children=tree['children']
        for child in children:
            children_of_root.append(child['name'])
        synonym_counts=create_mut_matrix(alphabet)
        non_synonym_counts=create_mut_matrix(alphabet)

    synonym_distribution=[0 for i in range(cds_start,cds_end)]
    non_synonym_distribution=[0 for i in range(cds_start,cds_end)]

    genetic_code=get_genetic_code()

    #counts the mutations, the root mutations are counted half in one direciton and half in the other
    with open(mutation_json_name) as augur:
        node_data = json.load(augur)
        nodes = node_data['nodes']
        for node in nodes:
            sequence=nodes[node]['sequence']
            mutations=nodes[node]['muts']
            for mutation in mutations:
                from_nt=mutation[0]
                to_nt=mutation[-1]
                position=int(mutation[1:-1])-1 #we remove one because our sequence starts from zero, while in the annotation it starts from 1.
                        
                if position>=cds_start and position<cds_end:
                    codon_offset=(position-cds_start)%3
                    to_codon=sequence[position-codon_offset:position+(3-codon_offset)] #use the sequence of the node to get the mutated sequence
                    if codon_offset==0:
                        from_codon=from_nt+to_codon[1:]
                    if codon_offset==1:
                        from_codon=to_codon[0]+from_nt+to_codon[2]
                    if codon_offset==2:
                        from_codon=to_codon[:2]+from_nt

                    if is_synonym(from_codon, to_codon, genetic_code):
                        synonym_distribution[position-cds_start]+=1
                        if node in children_of_root:
                            synonym_counts[from_nt+to_nt]+=1/2
                            synonym_counts[to_nt+from_nt]+=1/2
                        else:
                            synonym_counts[from_nt+to_nt]+=1
                    else:
                        non_synonym_distribution[position-cds_start]+=1
                        if node in children_of_root:
                            non_synonym_counts[from_nt+to_nt]+=1/2
                            non_synonym_counts[to_nt+from_nt]+=1/2
                        else:
                            non_synonym_counts[from_nt+to_nt]+=1
    
    total_synonym_mut=0
    for count in synonym_counts.values():
        total_synonym_mut+=count
    total_non_synonym_mut=0
    for count in non_synonym_counts.values():
        total_non_synonym_mut+=count

    return {"syn_distribution": synonym_distribution, "nonsyn_distribution":non_synonym_distribution, 
            "syn_counts":synonym_counts, "nonsyn_counts": non_synonym_counts,
            "tot_syn_mut":total_synonym_mut, "tot_nonsyn_mut":total_non_synonym_mut}

if __name__=="__main__":
    res = get_mutation_rates("wnv/results/ancestral_node_data.json")
    synonym_counts = res["syn_counts"]
    non_synonym_counts = res["nonsyn_counts"]
    print('synonym')
    print(synonym_counts)
    print('non synonym')
    print(non_synonym_counts)