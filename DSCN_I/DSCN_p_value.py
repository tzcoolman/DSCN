###############################
#required packages
###############################
import os,sys
import networkx as nx
from scipy.stats import norm
from scipy.stats import ttest_ind_from_stats
import time
###############################
'''
aaa=nx.Graph()
aaa.add_node('a',weight=5)
aaa.add_node('b',weight=9)
aaa.add_node('c',weight=4)
aaa.add_node('d',weight=3)
aaa['a']['b']['weight']=10
aaa['a']['c']['weight']=-3
aaa['c']['d']['weight']=5
#aaa.add_edge('a','b',weight=10)
#aaa.add_edge('b','c',weight=-3)
#aaa.add_edge('c','d',weight=5)
'''
#PAAD -0.015416104151 0.112180666564 7586924
def fit_ground_dis(graph,mu,std):
	total_regulation=[]
	#count=len(graph.edges())*2
	#print ("count:",count)
	if mu*std==0:
		for ele in graph.edges():
			#print ele[0],ele[1]
			#print graph.node[ele[0]]['weight']
			#print graph.node[ele[1]]['weight']
			edge_weight=graph[ele[0]][ele[1]]['weight']
			#print edge_weight
			node1_weight=graph.node[ele[0]]['weight']
			node2_weight=graph.node[ele[1]]['weight']
			#node1_weight=sum(graph.node[ele[0]]['weight'])/len(graph.node[ele[0]]['weight'])
			#node2_weight=sum(graph.node[ele[1]]['weight'])/len(graph.node[ele[1]]['weight'])
			total_regulation.append(node1_weight*edge_weight)
			total_regulation.append(node2_weight*edge_weight)
		mu, std = norm.fit(total_regulation)
	return (mu,std)

def calculate_p(sample_mu,sample_std,sample_num,pop_mu,pop_std,pop_num):
	return ttest_ind_from_stats(sample_mu,sample_std,sample_num,pop_mu,pop_std,pop_num)[1]
	#t_stat=result[0]
	#p=result[1]
	#return p

#print (calculate_p(5,20,10,70,49,2000))
if __name__ == '__main__':
	#print fit_ground_dis(aaa)
	print (calculate_p(-0.0364,0.22,200,-0.0154,0.112,7586924))
