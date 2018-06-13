### July 29, 2017
### DrugComboExplorer Version 1.0
### Author: Lei Huang
### Compute the synergistic drug combinations using pathway based similarity score

import csv
import numpy as np
import os
import math
import re
import sys
import time

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


drug_gene_reposition = {}

drug_gene_name = []

### First, we are going to load the repositioned signatures of drugs;

pathway_no = 242

res_folder = []

for i in range(1,2):

    res_folder.append("./oci_ly3_group_pathway_target_effect/")
    

drug_gene_name_total = []

drug_gene_dose_name_total = []

for i in range(1,pathway_no+1):
    drug_gene_reposition[i]=[]

for j in range(0,1):
    for i in range(1,pathway_no+1):

        oci_ly3_cell_line = open(res_folder[j] + 'oci_ly3_score_drug_refine_pathway_res_' + str(i) +'.txt','rU')
    
        oci_ly3_read = csv.reader(oci_ly3_cell_line, delimiter='\t')
    
        drug_name_profile_temp = []
        drug_gene_repo_temp = []
        drug_gene_name_temp = []
        drug_gene_dose_name_temp = []
        
        for item in oci_ly3_read:
            drug_gene_repo_temp.append(map(float,item[2:(len(item)-1)]))                                 
            drug_gene_dose_name_temp.append(item[1])
            drug_gene_name_temp.append(item[0])
            
        drug_gene_reposition[i] = drug_gene_reposition[i] + drug_gene_repo_temp
    drug_gene_name_total = drug_gene_name_total + drug_gene_name_temp
    drug_gene_dose_name_total = drug_gene_dose_name_total + drug_gene_dose_name_temp

for i in range(1,pathway_no+1):
    drug_gene_reposition[i]=np.array(drug_gene_reposition[i])
    
    #drug_gene_name[i] = drug_name_profile_temp

### This part we are going to compute the synergistic socre of the related drug pairs;

### we do not consider the same drug combiantions, so this kind of drug combinations should be removed.

drug_gene_name_unique = unique(drug_gene_name_total)

drug_gene_res = list(findInstances(drug_gene_name_unique,drug_gene_name_total))

drug_specific_group = {}

### we divide the same drug into specific drug groups, with drug in the same group, we have the same drug,
### drugs in the same group, we have the same drug name, but with different drug concentration.
### open the file of the pathway distances profiles

pathway_distance_f = open("./oci_ly3_pathway_distance_between_drugs_seed_genes_distance_test.txt","rU")

path_module_dist_a_b = {}
path_module_dist_b_a = {}

pathway_dist_read = csv.reader(pathway_distance_f)

for item in pathway_dist_read:
    item = item[0].split('\t')
    path_module_dist_a_b[item[0],item[1]] = float(item[2])
    path_module_dist_b_a[item[3],item[4]] = float(item[5])

for i in range(0,len(drug_gene_name_unique)):
    drug_specific_group[i] = drug_gene_res[i]

drug_combo_score_list = []

combination_score = []

for i in range(0,(len(drug_gene_name_unique)-1)):
    for j in range(i+1,len(drug_gene_name_unique)):
    #for j in range(0,len(drug_gene_name_unique)):
        for k in drug_specific_group[i]:
            #one_drug_dex = drug_specific_group[i][k]
            one_drug_dex = k
            one_drug_nodose_n = drug_gene_name_unique[i]
            one_drug_name = drug_gene_dose_name_total[k]
            for l in drug_specific_group[j]:
                #another_drug_dex = drug_specific_group[j][l]
                another_drug_dex = l
                another_drug_nodose_n = drug_gene_name_unique[j]
                another_drug_name = drug_gene_dose_name_total[l]
                sim_between_same_pathway_syner = 0.000
                sim_between_diff_pathway_syner = 0.000

                for rr in range(1,pathway_no+1):
                  
                    same_module_one = np.array(drug_gene_reposition[rr][one_drug_dex])
                    same_module_two = np.array(drug_gene_reposition[rr][another_drug_dex])
                    sim_module_score_numerator = same_module_one*same_module_two
                    
                    similarity = same_module_one * same_module_two
                    sim_value_mu_res = np.sum(abs(similarity))
                    
                    sim_value_mu_add = sim_module_score_numerator*sim_module_score_numerator
                    sim_value_mu = np.sum(sim_value_mu_add)
                    self_multiply_one = same_module_one * same_module_one
                    self_multiply_two = same_module_two * same_module_two
                    sim_module_score_a = np.sum(self_multiply_one)
                    sim_module_score_b = np.sum(self_multiply_two)
                    
                    
                    fen_mu_sim = sim_module_score_a + sim_module_score_b - sim_value_mu_res
                    if fen_mu_sim != 0.00:
                        sim_between_same_pathway = (sim_value_mu_res+0.00)*(sim_value_mu_res+0.00)/(sim_module_score_a + sim_module_score_b - sim_value_mu_res)
                        
                    else:
                        sim_between_same_pathway = 0.00 
                    
                    sim_between_same_pathway_syner = sim_between_same_pathway_syner + sim_between_same_pathway
                    

                diff_module_score = 0.000
                
                for tt in range(1,pathway_no):
                    for uu in range(tt+1,pathway_no+1):
                      
                        diff_module_one = np.array(drug_gene_reposition[tt][one_drug_dex])
                        diff_module_two = np.array(drug_gene_reposition[uu][another_drug_dex])
                        
                        diff_module_self_multi_one = diff_module_one * diff_module_one
                        diff_module_self_multi_two = diff_module_two * diff_module_two
                        
                        diff_self_multi_one = np.sum(diff_module_self_multi_one)
                        diff_self_multi_two = np.sum(diff_module_self_multi_two)
                        
                        pathway_distance_tt_uu = path_module_dist_a_b[str(tt),str(uu)]
                        pathway_distance_uu_tt = path_module_dist_b_a[str(uu),str(tt)]
                        
                        temp_multi_tt_uu = diff_self_multi_one * pathway_distance_tt_uu
                        temp_multi_uu_tt = diff_self_multi_two * pathway_distance_uu_tt
                        
                        diff_module_score = diff_module_score + temp_multi_tt_uu + temp_multi_uu_tt


                combo_score = sim_between_same_pathway_syner/(pathway_no) + diff_module_score*2/(pathway_no*(pathway_no-1))
                
                drug_combo_score_list.append([one_drug_nodose_n,another_drug_nodose_n,one_drug_name,another_drug_name,sim_between_same_pathway_syner,diff_module_score, combo_score])

                combination_score.append(combo_score)
  

specific_folder = time.strftime("%Y%m%d%d%d%d%d") + '_drug_combinations_oci_ly3_pathway_targeted_effects'

path = "./" + specific_folder

if not os.path.exists(path):
    os.makedirs(path)
    
drug_combination_res = open("./" + specific_folder + "/drug_combinations_oci_ly3_ranking"+".csv","wb")

drug_combo_res = csv.writer(drug_combination_res)

#### read the golden standard data from in vitro experiments

readgolden = open("./drug_synergy_data_IC20.txt","rU")

golden_res = csv.reader(readgolden, delimiter="\t")

combinations = []

goldenrank = []

goldenscore = []

rankindex = 0

removei = 0

for item in golden_res:
    if removei != 0:
        rankindex = rankindex + 1
        combinations.append(item[0].upper()+'#'+item[1].upper())
        goldenrank.append(rankindex)
        goldenscore.append(item[2])
    removei = removei+1
    
for item in drug_combo_score_list:
    drug_combo_res.writerow((item))

drug_combination_res.close()

unique_drug_combo_res = []

unique_drug_combo_score = []
    
for i in range(0,len(drug_combo_score_list)):
    unique_drug_combo_res.append(drug_combo_score_list[i][:2])
    unique_drug_combo_score.append(drug_combo_score_list[i][6])

unique_drug_combo_rank = unique(unique_drug_combo_res)

unique_drug_combo_score = np.array(unique_drug_combo_score)
    
refer_comb = list(findInstances(unique_drug_combo_rank,unique_drug_combo_res))

final_comb = []
final_average = []

for ii in range(0,len(refer_comb)):
    score = 0.0
    averagescore = 0.000
    for jj in range(0,len(refer_comb[ii])):
        #score = max(score,drug_combinations_alster[refer_comb[ii][jj]])
        averagescore = averagescore + float(unique_drug_combo_score[refer_comb[ii][jj]])
        score = max(score,float(unique_drug_combo_score[refer_comb[ii][jj]]))
    averagescore = averagescore/len(refer_comb[ii])
    final_comb.append(unique_drug_combo_rank[ii]+[score])
    final_average.append(unique_drug_combo_rank[ii]+[averagescore])

drugcombinationlist = []

combinationforrank = []

originaldrugindex = []

averageforrank = []

for i in range(0,len(final_comb)):
    drug_one = final_comb[i][0]
    drug_two = final_comb[i][1]
    if drug_one=='H-7, DIHYDROCHLORIDE':
        drug_one='H-7'
    elif drug_two=='H-7, DIHYDROCHLORIDE':
        drug_two='H-7'
    else:
        pass
        
    if drug_one == 'DOXORUBICIN HYDROCHLORIDE':
        drug_one = 'DOXORUBICIN'
    elif drug_two== 'DOXORUBICIN HYDROCHLORIDE':
        drug_two = 'DOXORUBICIN'
    else:
        pass
        
    if drug_one+'#'+drug_two in combinations:
        drugcombinationlist.append([drug_one,drug_two])
        combinationforrank.append(final_comb[i][2])
        averageforrank.append(final_average[i][2])
        originaldrugindex.append(combinations.index(drug_one+'#'+drug_two))
    elif drug_two+'#'+drug_one in combinations:
         drugcombinationlist.append([drug_two,drug_one])
         combinationforrank.append(final_comb[i][2])
         averageforrank.append(final_average[i][2])
         originaldrugindex.append(combinations.index(drug_two+'#'+drug_one))
    else:
         print 'key matching error',final_comb[i],combinations[i]

#comborank = sorted(range(len(combination_score)), key = lambda k:combination_score[k])
            
originalindex = np.array(originaldrugindex)

rankcombination = np.array(combinationforrank)

rankaveragecombination = np.array(averageforrank)

drugcombinationlists = np.array(drugcombinationlist)

indexsorting = originalindex.argsort()

drugcombinationoriginal = drugcombinationlists[indexsorting]

combinationrankings = rankcombination[indexsorting]

combinationaverageranking = rankaveragecombination[indexsorting]

finalcomborank = sorted(range(len(combinationrankings)), key = lambda k:combinationrankings[k],reverse=True)

finalavercomborank = sorted(range(len(combinationaverageranking)), key = lambda k:combinationaverageranking[k],reverse=True)

finalrank = map(lambda k:k+1, finalcomborank)

finalaverrank = map(lambda k:k+1, finalavercomborank)

arrayfinalrank = np.array(finalrank)

arrayfinalaverrank = np.array(finalaverrank)

### ranking the combination results

rankingcombinationres = rankcombination.argsort()[::-1]
rankedcombo = drugcombinationlists[rankingcombinationres]
rankedscores = rankcombination[rankingcombinationres]

rankingavercombinationres = rankaveragecombination.argsort()[::-1]
rankedavercombo = drugcombinationlists[rankingavercombinationres]
rankedaverscores = rankaveragecombination[rankingavercombinationres]


compare_drug_combinations_final_comb = open('./'+ specific_folder + '/drug_combinations_oci_ly3_ranking_compare_golden.csv','wb')

drug_combines_rank = csv.writer(compare_drug_combinations_final_comb)

for ii in range(0,len(drugcombinationoriginal)):
    #drug_combines_rank.writerow([drugcombinationoriginal[ii][0],drugcombinationoriginal[ii][1],finalrank.index])
    drug_combines_rank.writerow([drugcombinationoriginal[ii][0],drugcombinationoriginal[ii][1],list(np.where(arrayfinalrank==(ii+1))[0])[0]+1])

compare_drug_combinations_final_comb.close()

drug_combinations_final_comb = open('./'+ specific_folder + '/drug_combinations_for_oci_ly3_ranking_results.csv','wb')

drug_combines = csv.writer(drug_combinations_final_comb)

for jj in range(0,len(rankedcombo)):
    drug_combines.writerow([rankedcombo[jj][0],rankedcombo[jj][1],rankedscores[jj]])

drug_combinations_final_comb.close()

