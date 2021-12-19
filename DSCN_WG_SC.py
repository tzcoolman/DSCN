import os,sys
from sklearn import metrics as MS
from sklearn.cluster import SpectralClustering as SC
import numpy as np
import random
def load_exp(filename):
	count=0
	exp_hash={}
	gene_name={}
	f1=open(filename)
	for line in f1.readlines():
		temp=line.replace("\n","").replace("\r","").split("\t")
		name=temp[0]
		del temp[0]
		exp_hash[name]=temp
		gene_name[count]=name
		count+=1
	return (exp_hash,gene_name)

def load_fc(filename):
	fc_hash={}
	f1=open(filename)
	for line in f1.readlines():
		temp=line.replace("\n","").replace("\r","").split("\t")
		fc_hash[temp[0]]=float(temp[1])
	return fc_hash

def load_PPI(filename):
	PPI_hash={}
	f1=open(filename)
	for line in f1.readlines():
		temp=line[:-1].replace("\r","").split("\t")
		PPI_hash[temp[0]+"\t"+temp[1]]=1

exp_file="tumor_exp.txt.WG"
fc_file="tumor.FC.WG.nolog"
ppi_file="string_PPI.network.WG"
'''
exp_file="MAPK_PATHWAY/tumor_exp_overlap.FC.MAPK_pathway"
fc_file="MAPK_PATHWAY/tumor_exp_overlap.txt.MAPK_pathway"
ppi_file='MAPK_PATHWAY/PPI_overlap.network.MAPK_pathway'
'''
(exp_hash,gene_name)=load_exp(exp_file)
fc_hash=load_fc(fc_file)
#print len(exp_hash)
#print len(fc_hash)
#print len(gene_name)
ppi_hash=load_PPI(ppi_file)
cor_matrix=np.zeros((len(exp_hash),len(exp_hash)))
row_num=0
col_num=0
total_row=len(gene_name)
total_col=len(gene_name)
while row_num<total_row:
	col_num=0
	while col_num<total_col:
		if row_num==col_num:
			if row_num==col_num:
				cor_matrix[row_num][col_num]=fc_hash[gene_name[row_num]]
			elif cor_matrix[row_num][col_num]==0:
				if gene_name[row_num]+"\t"+gene_name[col_num] in ppi_hash or gene_name[col_num]+"\t"+gene_name[row_num] in ppi_hash:
					cor_matrix[row_num][col_num]=abs(np.corrcoef(exp_hash[gene_name[row_num]],exp_hash[gene_name[col_num]])[0,1])*(-1)
					cor_matrix[col_num][row_num]=cor_matrix[row_num][col_num]
		col_num+=1
	row_num+=1
row_num=0
col_num=0

while row_num<total_row:
	col_num=0
	row_sum=sum(cor_matrix[row_num])-cor_matrix[row_num][row_num]
        while col_num<total_col:
		if row_num!=col_num:
			if cor_matrix[row_num][col_num]!=0:
				cor_matrix[row_num][col_num]=cor_matrix[row_num][col_num]/row_sum*cor_matrix[row_num][row_num]
		col_num+=1
	row_num+=1
#SC_result=SC(random_state=10,n_clusters=800,gamma=10,assign_labels='discretize',affinity='precomputed',n_jobs=4,).fit_predict(cor_matrix)
np.random.seed(200)
#SC_result=SC(random_state=10,n_clusters=2,gamma=10,assign_labels='discretize',affinity='precomputed',n_jobs=4,).fit_predict(cor_matrix)
#SC_result2=SC(random_state=10,n_clusters=2,gamma=10,assign_labels='discretize',affinity='precomputed',n_jobs=4,).fit_predict(cor_matrix)
#print (MS.calinski_harabaz_score(cor_matrix, SC_result))
#print (MS.calinski_harabaz_score(cor_matrix, SC_result2))
#print len(SC_result.labels_)
#total=len(SC_result.labels_)
#count=0
#while count<total:
#	print (gene_name[count],SC_result.labels_[count])
#	count+=1
#X = np.array([[1, 1], [2, 1], [1, 0], [4, 7], [3, 5], [3, 6]])

def rearrange_cluster_result(cluster_array,gene_names):
	cluster_hash={}
	count=0
	total=len(cluster_array)
	while count<total:
		if cluster_array[count] in cluster_hash:
			cluster_hash[cluster_array[count]].append(str(gene_names[count]))
		else:
			cluster_hash[cluster_array[count]]=[str(gene_names[count])]
		count+=1
	return cluster_hash
if __name__=="__main__":
	max_score=0
	gamma_v=1.5
	SC_result=0
	max_sc_result=0
	while gamma_v<2:
		for seed in range(1,2):
			k=200
			total=300
			while k<total:
				try:
					SC_result=SC(n_init=k,random_state=seed,n_clusters=k,gamma=gamma_v,assign_labels='discretize',affinity='precomputed',n_jobs=8,).fit_predict(cor_matrix)
					CH_score=MS.calinski_harabaz_score(cor_matrix, SC_result)
					#print ("Seed", seed,"Calinski-Harabasz Score with gamma=", gamma_v, "n_clusters=", k,"score:", CH_score)
					#print (SC_result.labels_)
					if CH_score>max_score:
						max_score=CH_score
						max_sc_result=SC_result
						#SC(random_state=seed,n_clusters=k,gamma=gamma_v,assign_labels='discretize',affinity='precomputed',n_jobs=8,).fit(cor_matrix)
				except:
					print ("Error for K: ",k)
				k+=10
		gamma_v+=0.5
final_result=rearrange_cluster_result(max_sc_result,gene_name)
for key in final_result.keys():
	print str(key)+"\t"+"\t".join(final_result[key])
'''
	f3=open(sys.argv[1],'w')
	count=0
	total=len(max_sc_result)
	temp_x=""
	while count<total:
		temp_x+=gene_name[count]+"\t"+str(max_sc_result[count])+"\n"
		count+=1
	f3.write(temp_x[:-1])
	f3.close()
'''
#clustering = SC(n_clusters=2,n_components=40,assign_labels="kmeans",random_state=0).fit(X)
#print (clustering.labels_)
