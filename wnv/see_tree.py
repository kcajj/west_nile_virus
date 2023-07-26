from Bio import Phylo

#T_raw=Phylo.read('results/raw_tree.nwk','newick')
#T_raw.root_at_midpoint()
#Phylo.draw(T_raw)
T=Phylo.read('results/refined_tree.nwk','newick')
Phylo.draw(T)