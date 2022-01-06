###############################
#required packages
###############################
import os,sys
import numpy as np
import networkx as nx
import time
import multiprocessing as mp
###############################
#my code in other files
###############################
import DSCN_scorer_local_TN as DSCN_scorer
import DSCN_loading as DSCN_load
import DSCN_p_value as DSCN_p
def find_corr_cluster(clusters,target):
	for key in clusters.keys():
		if target in clusters[key]:
			return clusters[key]
	return []

def sub_sample_low_group(target,omics_data_hash,data_type):
	mean=sum(omics_data_hash[target])/len(omics_data_hash[target])
	count=0
	high_pos=[]
	low_pos=[]
	for exp_value in omics_data_hash[target]: # gene1: v1,v2,v3,...
		if exp_value>=mean:
			high_pos.append(count) #record sample position
		else:
			low_pos.append(count) #record sample position
		count+=1
	if data_type=="essen":
		return high_pos
	else:
		return low_pos
###################################################
#selected subgroups based on selected sample number
###################################################
def sample_by_num(sample_num_list,sample_list): 
	temp_list=[]
	count=0
	for ele in sample_list:
		if count in sample_num_list:
			temp_list.append(ele)
	count+=1
	return temp_list
####################################################
#calculate new pairwise correlation using subgroups
####################################################
def new_pairwise_corr(sample_list,exp_hash,graph):
	for edge in graph.edges():
		temp_list1=sample_by_num(sample_list,exp_hash[edge[0]])
		temp_list2=sample_by_num(sample_list,exp_hash[edge[1]])
		graph[edge[0]][edge[1]]['weight_new']=np.corrcoef(np.array(temp_list1),np.array(temp_list2))[0,1]
###################################################
#calculate new essentiality score using subgroups
###################################################
def new_essentiality(sample_list,essen_hash,graph):
	for node in graph.nodes():	
		temp_list=sample_by_num(sample_list,essen_hash[node])
		graph.node[node]['weight_new']=sum(temp_list)/len(temp_list) #get new average essentiality value for each node in the subnetwork
############################################################
#load exp, gene essen and PPI
############################################################
#target_file="../SL_synlethDB_pancrea.WG.list"
#exp_file="../cellline_exp.txt.WG"
#essen_file="../cellline.FC.WG"
#ppi_file="../string_PPI.network.WG"
#cluster_file="../WG.cluster"

target_file=sys.argv[1]
exp_file=sys.argv[2]
essen_file=sys.argv[3]
ppi_file=sys.argv[4]
cluster_file=sys.argv[5]


'''
target_file="../string_target_unique.list"
exp_file="R_clustering/cellline_exp_overlap.txt"
essen_file="R_clustering/cellline_exp_overlap.FC"
ppi_file="R_clustering/string_PPI_selected.network"
cluster_file="R_clustering/whole_genome_1_cluster.txt"
'''
exp_hash={}
cluster_list={}
essen_hash={}
cell_graph=nx.Graph()
target_list={}
#load_target(target_file,target_list)
#load_cluster(cluster_file,cluster_list)
DSCN_load.load_omic_data(exp_file,exp_hash)
DSCN_load.load_omic_data(essen_file,essen_hash)
#print (len(exp_hash))
#print (len(essen_hash))i
#print essen_hash.keys()
#target_list=load_target(target_file)
DSCN_load.load_common_gene_to_graph(ppi_file,cell_graph)
DSCN_load.load_target(target_file,target_list,cell_graph)
DSCN_load.load_cluster(cluster_file,cluster_list,cell_graph)
############################################################
#only calculate correlations occured within subnetwork
############################################################
DSCN_load.add_weight_to_sub_network(cell_graph,exp_hash,essen_hash)
bg_mu=-0.015416104151
bg_std=0.112180666564
bg_stat=DSCN_p.fit_ground_dis(cell_graph,bg_mu,bg_std)
(bg_mu,bg_std,count)=bg_stat
print bg_mu,bg_sigma,count
############################################################
#print (list(cell_graph.subgraph(['POLE2','HMMR']).edges()))
############################################################
#score the first target
############################################################
def execute_task(task_set):
	result=""
	each_target1=task_set[0]
	t_list=task_set[1]
############################################################
##initialize subnetwork1
############################################################
	SN_1=find_corr_cluster(cluster_list,each_target1)
	SG_1=cell_graph.subgraph(list(SN_1))
############################################################
##score subnetwork1
############################################################
	target1_task_edges=DSCN_scorer.find_influ_neighbors(SG_1,each_target1) # gather all tasks (edges) to be examined. e.g. <target,source,order>
	(temp_score,score_list)=DSCN_scorer.score_influence(each_target1,SG_1,target1_task_edges,1) # score influence for first target
############################################################
##resampling after target1 being KN
############################################################
	each_target1_sub_pos=sub_sample_low_group(each_target1,exp_hash,'exp') # select subgroups (exp)
	each_target1_sub_essen=sub_sample_low_group(each_target1,essen_hash,'essen') # select subgroups (essen)
############################################################
##update new node and edge weights in subnetwork1
############################################################
	new_pairwise_corr(each_target1_sub_pos,exp_hash,SG_1) #update edge weights
	new_essentiality(each_target1_sub_essen,essen_hash,SG_1) #update node weights
############################################################
##p-value 
############################################################
	#calculate_p(np.mean(score_list),np.std(score_list),op_num_1,bg_mu,bg_std,count)
############################################################
##record target1 score
############################################################
	result+="T1: "+each_target1+" : "+str(temp_score)+"\n"
############################################################
##score target2 given target1 being KN
############################################################
	for each_target2 in t_list:
		adj_type=1
############################################################
###intialize subnetwork2=(subnetwork1|target1)
############################################################
		SN_2=""
		SG_2=""
		if each_target2 != each_target1: #not the sample targets
			if each_target2 in SN_1:
				SG_2=SG_1
				adj_type=2
			else:
				SN_2=find_corr_cluster(cluster_list,each_target2)
				SG_2=cell_graph.subgraph(SN_2)
#############################################################
###score subnetwork2
#############################################################
			target2_task_edges=DSCN_scorer.find_influ_neighbors(SG_2,each_target2)
			(temp_score2,score_list2)=DSCN_scorer.score_influence(each_target2,SG_2,target2_task_edges,adj_type)
#############################################################
###p-value
#############################################################
			temp_total_list=score_list+score_list2
			calculate_p(np.mean(temp_total_list),np.std(temp_total_list),len(temp_total_list),bg_mu,bg_std,count)
#############################################################
###score 2nd target
#############################################################
			result+="T2: "+each_target1+" AND "+each_target2+": "+str(temp_score+temp_score2)
#############################################################
#return result
#############################################################
	return result[:-1]
#############################################################
#pushing first targets into list for parallelization
#############################################################
'''
if __name__ == '__main__':
        p=mp.Pool(8)
	total_task_set=[]
	count=0
	for target in target_list:
		total_task_set.append((target,target_list))
        total_result=p.map(execute_task,total_task_set)
	for each_chunk in total_result:
		print (total_result)
'''
f __name__ == '__main__':
        p=mp.Pool(8)
        total_task_set=[]
        count=0
        for target_combo in target_list:
                temp_combo=target_combo.split("\t")
                total_task_set.append((temp_combo[0],[temp_combo[1]]))
                total_task_set.append((temp_combo[1],[temp_combo[0]]))
        #for task in total_task_set:
        #       print task
        total_result=p.map(execute_task,total_task_set)
        for each_chunk in total_result:
                print (each_chunk)

