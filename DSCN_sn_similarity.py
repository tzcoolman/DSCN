import DSCN_loading as DSCN_load
import DSCN_scorer_local_TN as DSCN_scorer
import os,sys
import networkx as nx
#target_file="SL_synlethDB_pancrea.WG.list"
#t_exp_file="tumor_exp.txt.WG"
#t_fc_file="tumor.FC.WG"
#c_exp_file="cellline_exp.txt.WG"
#c_essen_file="cellline.FC.WG"
#ppi_file="string_PPI.network.WG"
#cluster_file="WG.cluster"

target_file=sys.argv[1]
t_exp_file=sys.argv[2]
t_fc_file=sys.argv[3]
c_exp_file=sys.argv[4]
c_essen_file=sys.argv[5]
ppi_file=sys.argv[6]
cluster_file=sys.argv[7]

'''
target_file="../string_target_unique.list"
exp_file="R_clustering/cellline_exp_overlap.txt"
essen_file="R_clustering/cellline_exp_overlap.FC"
ppi_file="R_clustering/string_PPI_selected.network"
cluster_file="R_clustering/whole_genome_1_cluster.txt"
'''
t_exp_hash={}
c_exp_hash={}
cluster_list={}
c_essen_hash={}
t_fc_hash={}
cell_graph=nx.Graph()
target_list={}

def make_tc_matrix(graph,t_exp_hash,t_fc_hash):
	row=len(graph.nodes())
	t_matrix=np.zeros(row,row)
	c_matrix=np.zeros(row,row)
	count_row=0
	for node1 in graph.nodes():
		count_col=0
		for node2 in graph.nodes():
			if count_row==count_col:
				c_matrix[count_row][count_row]=graph.nodes[node1]['weight']
				t_matrix[count_row][count_row]=t_fc_hash[node]
			else:
				if graph.has_edge(node1,node2):
					c_matrix[count_row][count_col]=graph[node1][node2]['weight']
					t_matrix[count_row][count_col]=abs(np.corrcoef(t_exp_hash[node1],t_exp_hash[node2])[0,1])
			count_col+=1
		count_row+=1
	return (t_matrix,c_matrix)

def normalize_matrix(func_matrix):
	abs_func_matrix=np.absolute(func_matrix)
	total_col,total_row=func_matrix.shape()
	row=0
	while row<total_row:
		col=0
		row_sum=sum(abs_func_matrix[row])-abs_func_matrix[row][row]
		while col<total_col:
			if row!=col:
				func_matrix[row][col]=func_matrix[row][col]/row_sum*abs_func_matrix[row][row]
			col+=1
		row+=1
	return func_matrix

def subn_similarity(cell_matrix,tissue_matrix):
	t_trace=tissue_matrix.trace()[0]
	c_trace=cell_matrix.trace()[0]
	total_col,total_row=cell_matrix.shape()
	row=0
	while row<total_row:
		cell_matrix[row][row]=cell_matrix[row][row]/c_trace*t_trace
		row+=1
	row=0
	cell_matrix=normalize_matrix(cell_matrix)
	tissue_matrix=normalize_matrix(tissue_matrix)
	difference=0.0
	while row<total_row:
		col=0
		while col<total_col:
			if col!=row:
				difference+=abs(cell_matrix[row][col]-tissue_matrix[row][col])
			col+=1
		row+=1
	return difference


if __name__ == '__main__':
######################################################
#load tissue data
######################################################
	tissue_graph=nx.Graph()
	DSCN_load.load_omic_data(t_exp_file,t_exp_hash)
	DSCN_load.load_omic_data(t_fc_file,t_fc_hash)
	DSCN_load.load_common_gene_to_graph(ppi_file,tissue_graph)
	load_cluster(cluster_file,cluster_list,tissue_graph)
	DSCN_load.add_weight_to_sub_network(tissue_graph,cluster_list)
######################################################
#load cellline data
######################################################
	cell_graph=nx.Graph()
	DSCN_load.load_omic_data(c_exp_file,c_exp_hash)
	DSCN_load.load_omic_data(c_essen_file,c_essen_hash)
	DSCN_load.load_common_gene_to_graph(ppi_file,cell_graph)
	#DSCN_load.add_weight_to_sub_network(cell_graph,cluster_list)
	load_cluster(cluster_file,cluster_list,cell_graph)
	DSCN_load.add_weight_to_sub_network(cell_graph,cluster_list)
############################################################
#only calculate correlations occured within subnetwork
############################################################
	for subn in cluster.keys():
		t_matrix=make_tc_matrix(tissue_graph.subgraph(cluster[subn]),t_exp_hash,t_fc_hash)
		c_matrix=make_tc_matrix(cell_graph.subgraph(cluster[subn]),c_exp_hash,c_essen_hash)
		nt_matrix=normalize_matrix(t_matrix)
		nc_matrix=normalize_matrix(c_matrix)
		print ("cluster:",subn,subn_similarity(nc_matrix,nt_matrix))
