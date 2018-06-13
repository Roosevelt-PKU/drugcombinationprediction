### July 29, 2017
### DrugComboExplorer Version 1.0
### Authoer : Lei Huang
### Compute the distance of the proteins or pathway proteins pairs

import csv
import numpy
import os
import math
import re
import sys
import scipy
import networkx as nx

inf_v = 100000

pathway_no = 242

path_num = 242

target_ppi_net = nx.Graph()

ppi_net = nx.Graph()

def unique(seq, keepstr=True):
  t = type(seq)
  if t==str:
    t = (list, ''.join)[bool(keepstr)]
  seen = []
  return t(c for c in seq if not (c in seen or seen.append(c)))

def findInstances(list1, list2):
    """For each item in list1,
    return a list of offsets to its occurences in list2
    """
    for i in list1:
        yield [pos for pos,j in enumerate(list2) if i==j]


#### import the interactions from the pathway data

protein_interaction = open("./human_cancer_pathway_expr_signal_hprd_unique.txt","rU")

protein_inter = csv.reader(protein_interaction,delimiter=" ")

ppi_inter = []

for item in protein_inter:
    ppi_inter.append((item[0],item[1],1.00))

ppi_inter = unique(ppi_inter)

                     
##### compute the shortest path of the drug targets on the specific drug modules

drug_gene_reposition = {}
drug_gene_name = []

### This part we are going to load the targets of drugs;

p_threshold = 0.0001

recommendation_matrix = {}

drug_target_relation = open("./oci_ly3_drug_special_target.txt","rU")

drug_targets = csv.reader(drug_target_relation,delimiter="\t")

spy = 0

cmap_drug_name = []

protein_name = []

drug_target_relation = []

for item in drug_targets:
    #print item
    if spy == 0:
        cmap_drug_name = cmap_drug_name + item
    if spy != 0:
        protein_name.append(item[0])
        drug_target_relation.append(map(float,item[1:]))
    spy = spy +1


drug_target_relation = numpy.array(drug_target_relation)

protein_name = numpy.array(protein_name)

for i in range(0,len(cmap_drug_name)):
    target_th = numpy.where(drug_target_relation[:,i]>= p_threshold)
    related_protein = protein_name[target_th]
    recommendation_matrix[cmap_drug_name[i]] = related_protein

pathway_module = {}
expressed_gene_p = []

#### read the expressed genes for oci_ly3

lincs_expressed_gene = open("./oci_ly3_gene_expressed.txt","rU")

cmap_expr = csv.reader(lincs_expressed_gene,delimiter="\t")

for item in cmap_expr:
    #expressed_gene_p.append(item[0])
    expressed_gene_p.append(''.join(item))

for i in range(1, pathway_no+1):
    pathway_module[i] = []

files_current = os.listdir("./ly3_signal_info_update_seed/")

for j in range(1,path_num+1):
    read_file = files_current[j-1]
    open_hsnc_sig = open("./ly3_signal_info_update_seed/"+read_file,"r")
    hsnc_path_sig = csv.reader(open_hsnc_sig,delimiter="\t")
    for item in hsnc_path_sig:
        for subitem in item:
            if subitem in expressed_gene_p:
                pathway_module[j].append(subitem)

### get the small sub-network from the available gene list

pathway_module_gene_list = []

for jj in range(1,pathway_no+1):
    pathway_module_gene_list = pathway_module_gene_list + pathway_module[jj]

target_pathway = []

for jjj in range(0,len(ppi_inter)):
    if ((ppi_inter[jjj][0] in pathway_module_gene_list)&(ppi_inter[jjj][1] in pathway_module_gene_list)):
        target_pathway.append((ppi_inter[jjj][0],ppi_inter[jjj][1],1.00))
           
target_ppi_net.add_weighted_edges_from(target_pathway)

ppi_net.add_weighted_edges_from(ppi_inter)


read_files = open("./oci_ly3_pathway_gene_collect_for_seed_genes.txt","rU")                            

sub_net = csv.reader(read_files,delimiter="\t")

arlist_gene = []

arlist = []

for item in sub_net:
    arlist_gene = arlist_gene + [''.join(item)]
    arlist = arlist + [''.join(item)]
   
short_path_length = {}

shortest_p_r = open("./shortest_path_protein_oci_ly3_seed_genes.txt","rU")

short_f = csv.reader(shortest_p_r,delimiter="\t")

index = 0

shortest_path_len = []

protein_list = []

for item in short_f:
    if index == 0:
        protein_list.append(item)
    if index >= 1:
        shortest_path_len.append(map(float,item[1:]))
    index = index + 1

shortest_path_len = numpy.array(shortest_path_len)

protein_list = protein_list[0]

protein_dic = {}

for i in range(0,len(protein_list)):
    protein_dic[protein_list[i]] = i


drug_combo_pathway_aa_bb = {}
drug_combo_pathway_bb_aa = {}

pathway_dist = open("./oci_ly3_pathway_distance_between_drugs_seed_genes_distance_test.txt","wb")

for ll in range(1,pathway_no):
    for kk in range(ll+1,pathway_no+1):
        
        one_pathway_genes = pathway_module[ll]
        another_pathway_genes = pathway_module[kk]

        fen_zi_a = 0.00
        
        for pp in range(0,len(one_pathway_genes)):
            
            min_shortest_dis_a = 0.000
            final_shortest_dis_a = 0.000
            
            spy_flag = 0
            temp_num = 0
            
            for jj in range(0,len(another_pathway_genes)):
                try:
                    if (one_pathway_genes[pp] == another_pathway_genes[jj]):
                        temp_short_a = 0.00
                    else:
                        temp_short_a = shortest_path_len[protein_dic[one_pathway_genes[pp]],protein_dic[another_pathway_genes[jj]]]

                    spy_flag = 1
                    temp_num = temp_num + 1.0000
                except:
                    temp_short_a = 100000
                
                if (spy_flag == 1)&(temp_short_a!=100000):
                    min_shortest_dis_a = min_shortest_dis_a + temp_short_a
                if (spy_flag == 0):
                    final_shortest_dis_a = final_shortest_dis_a + temp_short_a

            if (spy_flag==1):    
                
                fen_zi_a = fen_zi_a + math.exp(-min_shortest_dis_a/(temp_num))/temp_num + 0.000

            else:
                fen_zi_a = fen_zi_a + math.exp(-final_shortest_dis_a)/temp_num + 0.000

        
        fen_zi_b = 0.00

        for kkk in range(0,len(another_pathway_genes)):
            
            min_shortest_dis_b = 0.00
            final_shortest_dis_b = 0.00
            spy_f = 0
            temp_num_b = 0
                        
            for jjj in range(0,len(one_pathway_genes)):
                
                try:
                    if (another_pathway_genes[kkk] == one_pathway_genes[jjj]):
                        temp_short_b = 0.00
                    else:
                        temp_short_b = shortest_path_len[protein_dic[another_pathway_genes[kkk]],protein_dic[one_pathway_genes[jjj]]]
                        
                    spy_f = 1
                    temp_num_b = temp_num_b + 1.0000
                except:
                    temp_short_b = 100000

                if (spy_f == 1)&(temp_short_b!=100000):
                    
                    min_shortest_dis_b = min_shortest_dis_b + temp_short_b
                    
                if (spy_f == 0):
                    
                    final_shortest_dis_b = final_shortest_dis_b + temp_short_b

            if(spy_f==1):
                
                fen_zi_b = fen_zi_b + math.exp(-min_shortest_dis_b/(temp_num_b))/temp_num_b + 0.000
                                                                          
            else:
                fen_zi_b = fen_zi_b + math.exp(-final_shortest_dis_b)/temp_num_b + 0.000
                

        drug_combo_pathway_aa_bb[ll,kk,'a','b'] = fen_zi_a
        drug_combo_pathway_bb_aa[ll,kk,'b','a'] = fen_zi_b
        
        pathway_dist.write('%d\t%d\t'%(ll,kk))
        pathway_dist.write('%6.4f\t'%fen_zi_a)
        pathway_dist.write('%d\t%d\t'%(kk,ll))
        pathway_dist.write('%6.4f\t'%fen_zi_b)
        pathway_dist.write('\n')

pathway_dist.close()








