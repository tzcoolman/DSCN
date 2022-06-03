import networkx as nx
import numpy as np
import DSCN_scorer_local_TN as DSCN_scorer

def copy_list(source,target):
	target=[]
	for ele in source:
		target.append(ele)
	return target

def sibling_connection_check(nx_graph,sibling_set):
	result=[]
	for ele1 in sibling_set:
		for ele2 in sibling_set:
			if ele1!=ele2:
				if nx_graph.has_edge(ele1,ele2):
					result.append(ele1,ele2)
	return result

def single_weight(nx_graph,cluster_id,cluster,FC,d_targets,corr_matrix): #load edge weight
        cluster_list=cluster.split("\t")
        temp_graph=nx.subgraph(graph,cluster_list)
        for d_target in d_targets.keys():
                temp_total_dmg=FC[d_target]
                previous_node=""
                if d_target not in cluster_list:
                        continue
                for ele in cluster_list:
                        temp=1.0
                        #print d_target+":"+ele
                        #print '\t'.join(cluster_list)
                        #print temp_graph.edges()
                        s_path=[]
                        try:
                                #print d_target+":"+cluster_list
                                s_path=nx.dijkstra_path(temp_graph,d_target,ele)
                                saved_dijkstra_path[d_target+"\t"+ele]=s_path
                        except:
                                s_path=[]
                        if s_path==[]:
                                continue
                        else:
                                previous_node=d_target
                                for node in s_path:
                                        temp*=corr_matrix[previous_node+"\t"+node]
                                        previous_node=node
                        temp*=FC[previous_node]
                        rank_result[d_target]=temp_total_dmg
                        temp_total_dmg+=temp
                #print d_target+":"+str(temp_total_dmg)+":"+cluster_id
                #rank_result[d_target+"\t"+s_path[-1]]=temp_total_dmg
        #print dmg_score+"\t"+b_d_target
def find_order_neighbors(nx_graph,target_node):
        count=0
        neighbor_set=[]
        for ele in nx_graph.neighbors(target_node):
                neighbor_set.append((target_node,ele))
        return neighbor_set

def find_influ_neighbors(nx_graph,target_node):
	gene_set=[target_node]
	temp_graph=nx_graph.copy()
	count=0
	neighbors=[]
	while (count<1000000):
		source_set=[]
		#print gene_set
		for gene in gene_set:
			count=1
			temp=temp_graph.neighbors(gene)
			if temp==[]:
				#nx_graph.remove_node(gene)
				continue
			for ele in temp:
				count=1
				if ele not in gene_set:
					if ele not in source_set:
						source_set.append(ele)
				else:
					count=2
				neighbors.append([ele,gene,count])
				neighbors.append([gene,ele,count])
			temp_graph.remove_node(gene)
		#print "------------"
		#print neighbors
		#print "------------"
		count+=1
		if source_set==[]:
			break
		else:
			#sibling_connection_check(edge_list,gene_set)
			gene_set=copy_list(source_set,gene_set)
	return neighbors

def score_influence(target_node,graph,target_edge_list,order):
	itself_score=graph.node[target_node]['weight'] #average essentiality value
	'''
	if order==1:
		weight=nx.get_edge_attributes(graph,'weight') #original weight
	else:
		weight=nx.get_edge_attributes(graph,'weight_new') #adjusted weight
	#print target_edge_list
	#print weight[(each_task[1],each_task[0])]
	'''
	total_score=0
	target_node_weight=0
	edge_weight=0
	gene_target_weight=0
	gene_source_weight=0
	operation_num=0
	temp_var_list=[]
	if order==1:
		target_node_weight=graph.node[target_node]['weight']
	else:
		target_node_weight=graph.node[target_node]['weight_new']
	#print (target_node_weight)
	#for each_node_weight in target_node_weight:
	total_score+=np.mean(target_node_weight)
	operation_num+=1
	temp_var_list.append(target_node_weight)
	for each_task in target_edge_list:
		#print '0'
		if order==1:
			edge_weight=graph[each_task[0]][each_task[1]]['weight'] #original edge weight
			gene_target_weight=np.mean(graph.node[each_task[0]]['weight']) #original essentiality 
			#print (each_task[0])
			gene_source_weight=graph.node[each_task[1]]['weight']
		else:
			edge_weight=graph[each_task[0]][each_task[1]]['weight_new'] #adjusted weight from subsampling
			gene_target_weight=np.mean(graph.node[each_task[0]]['weight_new']) #adjusted essentiality
			gene_source_weight=np.mean(graph.node[each_task[1]]['weight_new'])
		#print (edge_weight)
		total_score+=float(gene_target_weight)*float(edge_weight)
		operation_num+=1
		temp_var_list.append(gene_target_weight*edge_weight)
		if each_task==2:
			total_score+=gene_source_weight*edge_weight
			operation_num+=1
			temp_var_list.append(gene_target_weight*edge_weight)
		#itself_score+=weight[(each_task[0],each_task[1])]*(sum(target_hash[each_task[0]])/len(target_hash[each_task[0]])) #edge weight times mean essentiality value
	return (total_score,temp_var_list)
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
#print find_influ_neighbors(dick,1)
