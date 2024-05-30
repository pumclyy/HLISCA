# cellchat

CellChat applies a signaling molecule interaction database (CellChatDB.human) to predict intercellular communication patterns based on differentially overexpressed receptors and ligands, including soluble agonists/antagonists as well as membrane-bound receptors/co-receptors. Specifically, a significant cell-cell communication was detected at the level of ligand–receptor pairs, based on a statistical test that randomly permutes the group labels of cells (permutation *P* < 0.05) and then recalculates the communication probability. Known ligand–receptor pairs from the CellChatDB database were used to compute the communication probability at the pair level and then summarized into communication strength at the pathway level. The [KEGG pathway database](https://www.genome.jp/kegg/pathway.html) was used to infer the relationship between ligand–receptor pairs and pathways.

First, we splited our data into ever-smokers and never-smokers to compute the communication probability independently. **(split.R)**

Then, we compute the communication probability and save the results of the communication probability. By comparing the probability between ever-smokers and never-smokers, We nominated the top 5 pathways showing the largest difference in the relative sum of communication strength of all cell types between ever- and never-smoker groups in each direction (top 5 increased and top 5 decreased in ever-smokers). **(cellchat.R)**

We performed same things in the HLCA dataset to validate the trend of MHC-I and MHC-II communication
