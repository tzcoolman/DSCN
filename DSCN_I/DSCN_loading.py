import networkx as nx
import numpy as np

def load_omic_data(filename,data_hash):
	f1=open(filename)
	count=0
	for line in f1.readlines():
		temp_line=line.replace("\n","").replace("\r","").split("\t")
		name=temp_line[0]
		del temp_line[0]
		data_hash[name]=map(float,temp_line)
	f1.close()

def load_cluster(filename,clusters,nx_graph): #clustering result
        f1=open(filename) #cluster file
        for line in f1.readlines():
                temp_list=[]
                temp=line.replace("\n","").replace("\r","")
                if temp=="":
                        continue
                else:
                        temp=temp.split("\t")
                        c_name=temp[0]
                        del temp[0]
                        for ele in temp:
                                if ele in nx_graph:
                                        temp_list.append(ele)
                        clusters[c_name]=temp_list
        f1.close()
'''
def load_target(filename,target_list,nx_graph):
        #target_list={}
        f1=open(filename)
        for line in f1.readlines():
                temp=line.replace("\n","").split("\t")
                if temp[0] in nx_graph:
                        target_list[temp[0]]=1
        f1.close()
'''
def load_target(filename,target_list): #load drug target
	f1=open(filename)
	for line in f1.readlines():
		temp=line.replace("\n","").replace("\r","").split("\t")
		target_list[temp[0]]=1
	f1.close()

def load_common_gene_to_graph(PPI_file,graph):
        f1=open(PPI_file)
        for line in f1.readlines():
                temp=line.replace("\n","").replace("\r","").split("\t")
                graph.add_edge(temp[0],temp[1])
                #graph.add_edge(temp[0],temp[1],weight=np.corrcoef(np.array(exp_hash[temp[0]]),np.array(exp_hash[temp[1]]))[0,1])
        f1.close()
#-0.015416104151 0.112180666564 7586924
def add_weight_to_sub_network(graph,exp_hash,essen_hash):
	for node in graph.nodes():
		graph.node[node]['weight']=np.mean(essen_hash[node])
	for edge in graph.edges():
		graph[edge[0]][edge[1]]['weight']=np.corrcoef(np.array(exp_hash[edge[0]]),np.array(exp_hash[edge[1]]))[0,1]
		#graph.add_edge(temp[0],temp[1],weight=np.corrcoef(np.array(exp_hash[edge[0]]),np.array(exp_hash[edge[1]]))[0,1])
        #for cluster in cluster_list.values():
        #        temp_graph=graph.subgraph(cluster)
        #        for edge in temp_graph.edges():
        #                graph[edge[0]][edge[1]]['weight']=np.corrcoef(np.array(exp_hash[edge[0]]),np.array(exp_hash[edge[1]]))[0,1]
