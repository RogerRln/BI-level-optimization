# -*- coding: utf-8 -*-
"""
@author: rrodrigu

This function parses the gene_reaction_rules of the enzymatically-catalyzed rxns, identifying the genes involved in each reaction.
Later, the function calculates the gene expression ratio that will constrain each enzymatically-catalyzed rxn.
The function returns a vector of 1258 gene expression ratios, one ratio corresponding to one enzymatically-catalyzed rxn.

Usage:
ratios = Calculate_Gexratios('simulated_condition_label' , 'baseline_condition_label',  gene_exp_matrix)
"""

import os
import cobra
import numpy as np
import re

### Set your working directory #####
os.chdir('/path/to/directory/containing/the/below/files')


model= cobra.io.read_sbml_model("gb-2009-10-6-r69-s4.xml")  # Bacillus iBsu1103 metabolic model
peg_ids= np.genfromtxt("peg_seedsubti.txt", delimiter='\n', dtype=None) # Gene IDs used by Bacillus metabolic model
bsu_ids= np.genfromtxt("bsu_seedsubti.txt", delimiter='\n', dtype=None) # Gene IDs used by the Gene expression matrix
gene_exp= np.genfromtxt("gene_expression.csv", delimiter=',', dtype=None)  # Gene expression matrix containing 11 different conditions (Columns) and 3852 genes (Rows)
#################################################################################



def Calculate_Gexratios (condition_simulated, baseline_name,  gene_exp):
          
        and_genes = []
        or_genes = []
        flag_and = 0
        flag_or = 0
        ratios = []

        condition_tested = np.nonzero(gene_exp[0] == condition_simulated)[0][0]
        condition_baseline = np.nonzero(gene_exp[0] == baseline_name)[0][0]
        
        lista_target_reactions= [str(i) for i in model.reactions]
        # we calculate the 'gex_ratios' from the gene_reaction_rule of eachn enzymatically-catalyzed rxn
        for key in lista_target_reactions:
                x= str(key)
                reaction_id = x
                gene_rule= model.reactions.get_by_id (reaction_id).gene_reaction_rule
                searchObj = re.findall( r'peg\.\d+|or|and|MANUAL|OPEN|PROBLEM|SPONTANEOUS|GROWMATCH', gene_rule, re.M|re.I)
                if searchObj:
                    if r'peg' in searchObj[0]:
                        if len(searchObj) > 1:  # Enter to this part of the function if the reaction is catalyzed by an enzymatic complex or isoenzymes
                            for j, i in enumerate(searchObj):
                                if r'peg' in i:
                                    index = list(peg_ids).index(i)
                                    bsu = bsu_ids[index]
                                    
                                    if (bsu in gene_exp[:,1]):    
                                        gene_row = np.nonzero(gene_exp[:,1] == bsu)[0][0]
                                        gex_ratio = 2**(float(gene_exp[gene_row, condition_tested])) / 2**(float(gene_exp[gene_row, condition_baseline]))
                                    
                                    else:
                                        gex_ratio = 1  
                                        
                                elif i == 'and':
                                    and_genes.append(gex_ratio)
                                    flag_and= 1
                                    flag_or= 0                                    
                                elif i == 'or':
                                    if and_genes:
                                        and_genes.append(gex_ratio)
                                        mean_and = np.mean(np.array(and_genes))
                                        or_genes.append(mean_and)
                                        and_genes = []
                                        flag_or= 1
                                        flag_and= 0
                                    else:
                                        or_genes.append(gex_ratio)
                                        flag_or= 1
                                        flag_and= 0                                        
                                if i == searchObj[-1] and j == len(searchObj)-1:  ## If element 'i' is the last element of the 'searchObj' vector and there was an 'OR' or 'AND' before
                                    if flag_and:                                  ## we save the gex_ratio on the AND|OR vector
                                        and_genes.append(gex_ratio)
                                        flag_and=0
                                    elif flag_or:
                                        or_genes.append(gex_ratio)
                                        flag_or=0
                                                                                              
                            if and_genes and not or_genes:
                                gex_ratio = np.mean(np.array(and_genes))
                                    
                            elif or_genes and and_genes:
                                mean_and = np.mean(np.array(and_genes))
                                or_genes.append(mean_and)
                                gex_ratio = sum(or_genes)
                                    
                            else:
                                gex_ratio = sum(or_genes) 
                                    
                            ratios.append(gex_ratio)
                            and_genes= []
                            or_genes =[]   
                                
                           
                        elif len(searchObj) == 1:    ## If reaction is catalyzed by one enzyme, then enter to this part of the function    
                            index = list(peg_ids).index(searchObj[0])
                            bsu = bsu_ids[index]                          
                            if bsu in gene_exp[:,1]:
                                gene_row = np.nonzero(gene_exp[:,1] == bsu)[0][0]
                                gex_ratio = 2**(float(gene_exp[gene_row, condition_tested])) / 2**(float(gene_exp[gene_row, condition_baseline]))
                            else:
                                gex_ratio= 1
                            
                            ratios.append(gex_ratio)
        return(ratios)                    

                  
# Media conditions labels
                  
#'Fru_1' = Fructose
# 'Glu_1' = Glucose
# 'Gly_1' = Glycerol
# 'Glucon_1'= Gluconate
# 'Mal_1'  = Malate
# 'Pyr_1' = pyruvate
# 'M.G_1' = Malate + Glucose
# 'LBexp_1' = LB
# 'SMMPr_1' = SMM
# 'S0_1' = CH       
# 'G.S_1 = Glutamate + Succinate    
                  
media_conditions = ['Fru_1', 'Glu_1', 'Gly_1', 'Glucon_1', 'Mal_1', 'Pyr_1', 'M.G_1', 'LBexp_1', 'SMMPr_1', 'S0_1' , 'G.S_1']

# Example: Calculate Gene Expression Ratios for the simulated condition 'CH' using 'SMM' as baseline condition
#Usage:
#ratios = Calculate_Gexratios('simulated_condition_label' , 'baseline_condition_label',  gene_exp_matrix)

ratios = Calculate_Gexratios('S0_1', 'SMMPr_1', gene_exp )

# Vector containing the 1258 enzymatically-catalyzed rxn IDs
enz_rxns = [str(rxn) for rxn in model.reactions if r'peg' in model.reactions.get_by_id(str(rxn)).gene_reaction_rule]
