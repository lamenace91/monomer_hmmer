#! /usr/bin/env python

import sys
import re
import os
import subprocess
import csv
import numpy as np
from numpy import genfromtxt
from Bio import AlignIO
import time
import argparse
import tempfile
from Bio import SeqIO
import math
#from scipy import DataFrame
import pandas
from pandas import *
from multiprocessing import Lock, Process, Queue, current_process
import operator
from Bio.SeqRecord import SeqRecord
import my_function as mf


def get_rolling_alignments(seq_seq, mono_seq, beg, verbose=0):
	#############################
	rollingCircle=[]
	mm=beg-1
	ss=0
	for ii in range(0, len(mono_seq)):
		if mono_seq[ii] != "-" and seq_seq[ii] != "-" :
			mm=mm+1	
			ss=ss+1
			rollingCircle.append(mm)
		elif mono_seq[ii] == "-" and seq_seq[ii] != "-" :
			ss=ss+1
			rollingCircle.append(0)
		elif mono_seq[ii] != "-" and seq_seq[ii] == "-" :
			mm=mm+1			
	rollingCircle=np.array(rollingCircle)
	rollingCircle=rollingCircle
	return(rollingCircle)



######################################
## search best hit by using glsearch program (fasta package)
######################################
def	newhit_by_glssearch(query, library, beg=-1, end=-1, algo='gl'):
	options=" -m B -b 1 -d 1 -E 10000 " 
	region=""
	if beg > 0 and end > 0:
		region=":"+str(beg)+"-"+str(end)
	if beg > 0 and end <= 0:
		region=":"+str(beg)
		
	if algo == 'gl':
		command="glsearch36 "+options+" "+query+" "+library+region
	else:
		command="ssearch36 "+options+" "+query+" "+library+region
	hit={}
	hit['query_beg']=0
	hit['query_end']=0
	hit['subject_beg']=0
	hit['subject_end']=0
	hit['subject_strand']=0
	hit['ali_simil']=0.0
	hit['ali_length']=0
	hit['ali_score']=0.0
	hit['ali_eval']=0.0
	hit['query']=""
	hit['subject']=""
	result=""
	try:
		#print("------------------------>")
		result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
	except subprocess.CalledProcessError as exc:
		print 'error: code={}, out="{}"'.format(
        exc.returncode, exc.output,
        )
	lines=result.split('\n')
	nb_hit=0
	for line in result.split('\n'):
		words=line.split()
		if len(words) > 2 and words[0] == "Score" :
			nb_hit=nb_hit+1
		if len(words) > 7 and words[0] == "Score" and nb_hit == 1:
			hit['ali_score']=float(words[2])
			hit['ali_eval']=float(words[7])
		if len(words) > 4 and words[0] == "Identities" and nb_hit == 1:
			hit['ali_simil']=float(re.sub('[(),%]','',words[3]))
			hit['ali_length']=int(re.sub('.*/','',words[2]))	
		if len(words) > 2 and words[0] == "Query" and nb_hit == 1:
			if hit['query'] == "":
				hit['query_beg']=int(words[1])
			hit['query_end']=int(words[3])
			hit['query']=hit['query']+words[2]
		if len(words) > 2 and words[0] == "Sbjct" and nb_hit == 1:
			if hit['subject'] == "":
				hit['subject_beg']=int(words[1])
			hit['subject_end']=int(words[3])
			hit['subject']=hit['subject']+words[2]
	return(hit)
######################################
## Run analysis for one sequence
######################################
#def sequence_analysis(work_queue, done_queue, mono_file, minimal_length=10, minimal_simil=0.0, minimal_score=0.0, maximal_eval=10000000.0, verbose=0):
def sequence_analysis(file_name, mono_file, output_file_monomers, output_file_tabular, output_file_rc, rc_index=1, minimal_length=10, minimal_simil=0.0, minimal_score=0.0, maximal_eval=10000000.0, verbose=0):
    #print "Starting analysis process "+current_process().name+" ..."
    #for file_name in iter(work_queue.get, 'STOP'):
	#print "#########################################"
	#print current_process().name+" : "+file_name
        #sys.stdout.flush()
	t1=time.time()
	t1cpu=time.clock()
    
    
	seq_len=mf.get_seq_length(file_name)
	seq_name=seq_len.keys()[0]
	hits=[]
	hit={}
	hit['query_beg']=0
	hit['query_end']=0
	hit['subject_beg']=0
	hit['subject_end']=0
	hit['subject_strand']=0
	hit['ali_simil']=0.0
	hit['ali_length']=0
	hit['ali_score']=0.0
	hit['ali_eval']=0.0
	hit['query']=""
	hit['subject']=""
	hits.append(hit)
	hit={}
	hit['query_beg']=0
	hit['query_end']=0
	hit['subject_beg']=0
	hit['subject_end']=0
	hit['subject_strand']=0
	hit['ali_simil']=0.0
	hit['ali_length']=0
	hit['ali_score']=0.0
	hit['ali_eval']=0.0
	hit['subject_beg']=seq_len[seq_name]+1
	hit['query']=""
	hit['subject']=""
	hits.append(hit)
		
	ii=1
	#print seq_name
	while ii < len(hits):
		#print "##############"
		#for hh in hits:
			#print("%d %d - %d %d %d - %d %f" % (hh['query_beg'], hh['query_end'] ,hh['subject_beg'], hh['subject_end'], hh['subject_strand'], hh['ali_length'], hh['ali_score']))
		##print ii
		#print len(hits)
		#print seq_len
		if (hits[ii]['subject_beg']-hits[ii-1]['subject_end']) > 10:
			nh=newhit_by_glssearch(mono_file, file_name, hits[ii-1]['subject_end']+1, hits[ii]['subject_beg']-1, algo="gl")
			#print "--------------"
			#print hits[ii-1]['subject_end']+1
			#print hits[ii]['subject_beg']-1
			#print nh
			if nh['ali_simil']== 0.0:
				nh=newhit_by_glssearch(mono_file, file_name, hits[ii-1]['subject_end']+1, hits[ii]['subject_beg']-1, algo="s")
				#print "short ...."
				#print nh
			if nh['ali_simil'] >= minimal_simil and nh['ali_score'] >= minimal_score and nh['ali_length'] >= minimal_length and nh['ali_eval'] <= maximal_eval:
				hits.append(nh)
				hits.sort(key=operator.itemgetter('subject_beg'))
				hits=sorted(hits, key=lambda k: int(k['subject_beg'])) 
				ii=ii-1	
		ii=ii+1		
		
#	print "finished "	
	if verbose > 0:
		for hh in hits:
			print("%d %d - %d %d %d - %d %f" % (hh['query_beg'], hh['query_end'] ,hh['subject_beg'], hh['subject_end'], hh['subject_strand'], hh['ali_length'], hh['ali_score']))
	rc=np.array([])
	rc=np.array([])
	tmp=1
	for ii in range(1,len(hits)-1):
		for jj in range(tmp, hits[ii]['subject_beg']):
			rc=np.concatenate((rc, [-1]))
		gr= get_rolling_alignments(hits[ii]['subject'],hits[ii]['query'], hits[ii]['query_beg'], verbose)
		rc=np.concatenate((rc, gr))
		tmp=hits[ii]['subject_end']+1
	for jj in range(tmp, seq_len[seq_name]+1):
		rc=np.concatenate((rc, [-2]))
	
	writer_to_file(file_name, seq_name, hits, rc, output_file_monomers, output_file_tabular, output_file_rc, rc_index)
	#done_queue.put((file_name, seq_name, hits, rc))
	t2=time.time()
	t2cpu=time.clock()
	print "Ending analysis process: "+current_process().name+" file name: "+file_name+" seq name: "+seq_name+ " seq length:  "+str(seq_len[seq_name])+" wall-clock time: "+str(t2-t1)+" cpu time: " +str(t2cpu-t1cpu) + " !!!"
	return True























##########################################
#def writer_to_file(done_queue,output_file_monomers, output_file_tabular, output_file_rc, rc_index=1, verbose=0, naming=1):
def writer_to_file(seq_file, seq, m, rc, output_file_monomers, output_file_tabular, output_file_rc, rc_index=1, verbose=0, naming=1):
   #print "Starting writing process ..."
    #while 1:
        #seq_file, seq, m, rc = done_queue.get()
	#if seq_file != 'STOP':
		handle=open(output_file_tabular, 'a')
		print_analysis(seq, m)
		if verbose > 5:
			print ("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
			print(m)
			print ("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
		print_ascii(seq,m)
		for hh in m:
			if hh['query_beg'] != 0 or  hh['query_end'] != 0 :
				handle.write("%s - %d %d - %d %d %d - %d %f %f - %s %s\n" % (seq, hh['query_beg'], hh['query_end'] ,hh['subject_beg'], hh['subject_end'], hh['subject_strand'], hh['ali_length'], hh['ali_score'], hh['ali_eval'], hh['subject'], hh['query']))
		handle.close()
		
		# wrinting rollong_circle information
		if output_file_rc != "":
			handle=open(output_file_rc, 'a')
			handle.write(seq +" ")
			for cc in rc:
				handle.write( "," + str(int(cc)))
			handle.write("\n")
			handle.close()
			
		ss=mf.get_sequence(seq,seq_file).seq
		if 	len(ss) != len(rc):
				print "Error: "+seq+ " "+str(len(ss))+" "+str(len(rc))
				
		if output_file_monomers != "":		
			handle=open(output_file_monomers+str(rc_index), 'a')
			tmp=0
			ok=0
			yy=0
			for ii in range(0,len(ss)):
					yy=yy+1
					if ok == 0 or rc[ii] > 0 and (rc[ii] == rc_index or (tmp < rc_index and rc[ii] >=rc_index) or (tmp > rc[ii] and rc[ii] >=rc_index)):
						if ok == 1:
							handle.write( "\n")
#							print ""
##						print ">"+seq+"_"+str(ii)
						handle.write( ">"+seq+"_"+str(ii)+"\n")
						ok=1
					handle.write((ss[ii]))
##					print ss[ii],
					if rc[ii] > 0:
						tmp=rc[ii]
			handle.write("\n")
	##		print ""
			handle.close()
				
			os.remove(seq_file)
	#else:
		#break
		print "Ending writing process!!!"
		return True
#####################################
def split_file(seq_file, n=1):	
		split_file_list=[]
		if not os.path.exists(seq_file):
			return split_file_list
			
		handle=open(seq_file, "rU")
		nn=n
		first=True
		for record in SeqIO.parse(handle, "fasta") :
			if nn >= n:
				if not first:
					output_handle.close()
				seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
				output_handle=open(seq_tmpfile, "w")
				split_file_list.append(seq_tmpfile)
				first=False
			SeqIO.write(record, output_handle, "fasta")
			nn=nn+1
		return(split_file_list)

##########################################
def print_ascii(seq, m): 
	init_gap= m[1]['subject_beg']

	output=""
	if init_gap > 0:
		output=output+ ".."+str(init_gap)+".."
	for ii in range(1,len(m)-1):
		leng=m[ii]['subject_end']-m[ii]['subject_beg']
		if m[ii]['subject_strand'] == 1:
			output=output+ "--"+str( leng)+"-->" 
		else:
			output=output+ "<-"+str( leng)+"--" 
		dist=m[ii+1]['subject_beg']-m[ii]['subject_end']	
		if dist != 1 :
			output=output+".."+str(dist)+".."
	print output
			
	for ii in range(1,len(m)-1):
			print ' ',	
######################################
def print_analysis(seq, data):
	total_gap=0
	#dt = np.dtype([np.str_, np.float, np.int, np.int, np.int, np.int, np.int, np.str_, np.float, np.int, np.int])
	
	#dt= np.dtype([('monomer', 'a100'), ('score', 'f8'), ('l', 'i4'), ('b', 'i4'), ('c', 'i4'), ('d', 'i4'), ('e', 'i4'),('strand', 'a100'),('h', 'f8')  ,('f', 'i4'), ('g', 'i4') ])
	#xx=np.asarray(data)
	#data=np.atleast_2d(xx)
	#data.dtype=dt
	nb_lines=len(data)-1
	seq_len=int(data[nb_lines]['subject_beg'])
	init_gap= data[1]['subject_beg']
	term_gap= seq_len-int(data[nb_lines-1]['subject_end'])
	
	for ii in range(1,nb_lines):
		total_gap+=(data[ii]['subject_beg']-data[ii-1]['subject_end']-1)
	total_gap+=term_gap
	print(seq+" seq len: "+str(seq_len)+" initial_gap: "+str(init_gap)+" term gap: "+str(term_gap)+" total gap: "+str(total_gap))
	return(data)


######################################
parser = argparse.ArgumentParser(description='Search for tandemly repeated monomers ...', epilog='This program searchs for monomer limits into long sequences (-s option) by using the Smith and Waterman algorithm as implemented in _water_ (package EMBOSS, needed) and the monomers sequences (-m option)')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-m', dest="mono_file", required=True)
parser.add_argument('-o', dest="output_file", required=True)
parser.add_argument('-O', dest="output_file_tab", required=True)
parser.add_argument('-r', dest="output_file_rc", default="")
parser.add_argument('-i', dest="rc_index", type=int, default=1)
parser.add_argument('-f', dest="number_of_seq", type=int, default=1)

parser.add_argument('-n', dest="naming",     type=int,         default=2)
#parser.add_argument('-proc',  dest="processors",       type=int,         default=4)
parser.add_argument('-verb',  dest="verbose",          type=int,              default=10)
parser.add_argument('-wml',  dest="waterman_minimal_length", type=int, default=10)
parser.add_argument('-wms',  dest="waterman_minimal_score",  type=float, default=100)
parser.add_argument('-wmS',  dest="waterman_minimal_simil",  type=float, default=0.0)
parser.add_argument('-wge',  dest="waterman_extension_gap_penality", type=float, default=1)
parser.add_argument('-wgo',  dest="waterman_opening_gap_penality",   type=float, default=10)




args = parser.parse_args()
if not os.path.exists(args.seq_file):
	print("Error: sequence file ("+args.seq_file+") not found !!")
	quit()
if not os.path.exists(args.mono_file):
	print("Error: sequence file ("+args.mono_file+") not found !!")
	quit()
if args.number_of_seq < 1:
	args.number_of_seq=1
	
if os.path.exists(args.output_file):
    os.remove(args.output_file)
if os.path.exists(args.output_file_tab):
    os.remove(args.output_file_tab)
if args.output_file_rc != "" and os.path.exists(args.output_file_rc):
    os.remove(args.output_file_rc)



if args.verbose > 3:
	print ("spliting sequence file ...")
files=split_file(args.seq_file, args.number_of_seq)
nb_files=len(files)
if args.verbose > 3:
	print ("Number of files: "+str(nb_files))


#workers = args.processors
#work_queue = Queue()
#done_queue = Queue()
#processes = []
#for ff in files:
	#work_queue.put(ff)	
	
#q = Process(target=writer_to_file, args=(done_queue,args.output_file,args.output_file_tab, args.output_file_rc, args.rc_index, args.verbose, args.naming))
#q.start()



for ff in files:
		sequence_analysis(ff, args.mono_file,args.output_file,args.output_file_tab, args.output_file_rc, args.rc_index, args.verbose)
		
#for w in xrange(workers):
        #p = Process(target=sequence_analysis, args=(work_queue, done_queue, args.mono_file,args.verbose))						
        #p.start()
        #processes.append(p)
        #work_queue.put('STOP')

#for p in processes:
        #p.join()


#done_queue.put(('STOP','STOP','STOP','STOP'))

#print "Waiting for the writing  process ..."
#q.join()
	
print "Bye bye !!!"

