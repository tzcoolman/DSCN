import os,sys
#sys.argv[-1]=sys.argv[-1].replace("\r","")
##########################################################
#your tumor exp and fc files
##########################################################
#tumor_exp=sys.argv[1]
#tumor_fc=sys.argv[2]
##########################################################
#Your cell-line exp and fc files
##########################################################
#cell_exp=sys.argv[3]
#cell_fc=sys.argv[4]
##########################################################
#drug bank target list 
##########################################################
#target_list=sys.argv[5]
##########################################################
#STRING PPI network 
##########################################################
#ppi_net=sys.argv[6]
##########################################################
#add Your thread ID 
##########################################################
#file_path = tumor_exp.replace(tumor_exp.split("/")[-1],"")
##########################################################
prefix="Overlap_"
final_overlap={}

def load_hash(filename,f_type):
	x_final={}
	f1=open(filename)
	for line in f1.readlines():
		temp=line.replace("\n","").replace("\r","").split("\t")
		name=temp[0]
		del temp[0]
		if f_type==1:
			x_final[name]="\t".join(temp)
		else:
			x_final[name]=1
			x_final[temp[0]]=1
	return x_final

def sub_select(filename,x_type,file_path):
	#print (file_path)
	#print ("here")
	#print (prefix)
	f1=open(filename)
	sub_string=""
	count=0
	for line in f1.readlines():
		temp=line.replace("\n","").split("\t")
		if x_type==1:
			if temp[0] in final_overlap:
				sub_string+=line.replace("\r","")
		else:
			if temp[0] in final_overlap and temp[1] in final_overlap:
				sub_string+=line.replace("\r","")
		count+=1
	f1.close()
	filename=filename.split("/")[-1]
	f2=open(file_path+prefix+filename,'w')
	f2.write(sub_string[:-1])
	f2.close()

def batch_overlap(tumor_exp,tumor_fc,cell_exp,cell_fc,target_list,ppi_net,file_path):
	#file_path=tumor_exp.replace(tumor_exp.split("/")[-1],"")
	file_path=file_path.replace("\r","").replace("\n","")
	if file_path[-1]!="/":
		file_path+="/"
	t_exp=load_hash(tumor_exp,1)
	t_fc=load_hash(tumor_fc,1)
	c_exp=load_hash(cell_exp,1)
	c_fc=load_hash(cell_fc,1)
	ppi=load_hash(ppi_net,2)
	for key in t_exp.keys():
		if key in t_fc:
			if float(t_fc[key])>0 and key in c_exp and key in c_fc and key in ppi:
				final_overlap[key]=1
	sub_select(tumor_exp,1,file_path)
	sub_select(tumor_fc,1,file_path)
	sub_select(cell_exp,1,file_path)
	sub_select(cell_fc,1,file_path)
	sub_select(ppi_net,2,file_path)
	sub_select(target_list,1,file_path)
if __name__=="__main__":
	batch_overlap(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],'./')