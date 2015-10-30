#!/usr/local/bin/env python3.4

'''
you need BLAST+ on your Linux
'''

query_fasta_file='new_fasta.fasta'
blast_type='blastp'

import os
import re
import multiprocessing
import threading

def run_blast():
	if 'blastp'==blast_type.lower():
		os.system('makeblastdb -in myva.txt -dbtype prot -parse_seqids -out myva.txt.db')
		os.system('blastp -query '+query_fasta_file+' -out '+query_fasta_file+'.blastp -db myva.txt.db -outfmt 6 -evalue 1e-5')
		os.system('rm myva.txt.db.*')
	elif 'blastx'==blast_type.lower():
		os.system('makeblastdb -in myva.txt -dbtype prot -parse_seqids -out myva.txt.db')
		os.system('blastx -query '+query_fasta_file+' -out '+query_fasta_file+'.blastx -db myva.txt.db -outfmt 6 -evalue 1e-5')
		os.system('rm myva.txt.db.*')

def get_whog():
	global whog_name_func_dict;whog_name_func_dict={}
	global whog_name_gene_dict;whog_name_gene_dict={}
	global whog_name_COG_dict;whog_name_COG_dict={}
	with open('whog.txt') as whog:
		for i in whog:
			if '[' in i:
				m=i.rstrip()
				COG_num=re.findall('\] (COG[0-9]*?) ',m)[0]
			elif ':' in i:
				n=i.split(':')[0].lstrip()
				for j in [x for x in i.split(':')[1].rstrip().split(' ') if x]:
					whog_name_func_dict[j]=m
					whog_name_gene_dict[j]=n
					whog_name_COG_dict[j]=COG_num

def get_fun():
	global fun_dict;fun_dict={}
	global class_dict;class_dict={}
	with open('fun.txt') as fun:
		for i in fun:
			if '[' in i:
				lin=i.lstrip().rstrip().split(' ',1)
				fun_dict[lin[0][1]]=lin[1]
				class_dict[lin[0][1]]={}

def deal_blast_result():
	global blast_dict;blast_dict={}
	a=open(query_fasta_file+'.'+blast_type,'a')
	a.write('20	20	20	')
	a.close()
	with open(query_fasta_file+'.'+blast_type) as blast_result:
		start=1
		for i in blast_result:
			if i.rstrip():
				m=i.split('	')[0];n=i.split('	')[1];o=float(i.split('	')[2])
				if float(o) >= 30.0:
					if start == 1:
						mm=m;nn=n;oo=o
						start+=1
					elif start != 1:
						if m == mm and o > oo:
							mm=m;nn=n;oo=o
						elif m != mm:
							blast_dict[mm]=nn
							mm=m;nn=n;oo=o

def generate_COGnum_result():
	COGnum_dict={}
	fout_COGnum=open(query_fasta_file+'.COGnum','w')
	for i,j in blast_dict.items():
		if j in whog_name_func_dict.keys():
			if whog_name_COG_dict[j] in COGnum_dict.keys():
				COGnum_dict[whog_name_COG_dict[j]].append(i)
			else:
				COGnum_dict[whog_name_COG_dict[j]]=[]
				COGnum_dict[whog_name_COG_dict[j]].append(i)
	for i,j in COGnum_dict.items():
		print(i)
		fout_COGnum.write(i+'	'+str(len(j))+'\n	')
		for k in j:
			fout_COGnum.write(k+'\n	')
		fout_COGnum.write('\n')
	fout_COGnum.close()

def generate_COGclass_result():
	fout_COGclass=open(query_fasta_file+'.COGclass','w')
	for i,j in blast_dict.items():
		if j in whog_name_func_dict.keys():
			for k in re.findall('\[([A-Z]*?)\]',whog_name_func_dict[j])[0]:
				class_dict[k][i]=j
	for i in class_dict:
		fout_COGclass.write(i+'	'+fun_dict[i]+'\n	')
		for j,k in class_dict[i].items():
			fout_COGclass.write(j+'	'+whog_name_func_dict[k]+'\n	')
		fout_COGclass.write('\n')
	fout_COGclass.close()

if __name__=='__main__':
	p = multiprocessing.Process(target=run_blast, args=())
	p.start()
	threads_list=[]
	threads_list.append(threading.Thread(target=get_whog,args=()))
	threads_list.append(threading.Thread(target=get_fun,args=()))
	for t in threads_list:
		t.setDaemon(True)
		t.start()
	p.join()
	deal_blast_result()
	q = multiprocessing.Pool()
	q.apply_async(generate_COGnum_result, args=())
	q.apply_async(generate_COGclass_result, args=())
	q.close()
	q.join()