#! /usr/bin/env python

import re
import timeit
import os
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
from Bio.SeqRecord import SeqRecord
from my_function import extract_sequence_to_newFile

##########################################################
def extract_sequence_to_newFile(seq, from_file, to_file, beg=None, end=None, option='a', name=0 
	#, reverseComplement=False
	):
	found=0
	record_index = SeqIO.index(from_file, "fasta")
	if seq in record_index:
		if beg is None:
			beg=1
		if end is None:
			end=len(record_index[seq].seq)
		output_handle = open(to_file, option)
		if name ==0:
			seqid=record_index[seq].id+":"+str(beg)+"-"+str(end)
		if name == 1:
			seqid=record_index[seq].id
		
		SeqIO.write(SeqRecord(record_index[seq].seq[beg:end],id=seqid), output_handle, "fasta")
		output_handle.close()		
		return 1
	else:
		return 0	
	return found
##########################################################

def blast(seq_file, mono_file, proc=1, word_size=4, output_tmpblastfile="", rmout=0):	
	# blasting ...
	run=0
	if output_tmpblastfile == "":
		run=1
		output_tmpblastfile=tempfile.NamedTemporaryFile(delete=False).name
		
	output_tmplog=tempfile.NamedTemporaryFile(delete=False).name
	if not os.path.exists(seq_file):
    		return -1
	if not os.path.exists(mono_file):
    		return -2
	if rmout == 0:	
		print "files are not deleted ..."
		print "blast output: "+ output_tmpblastfile
		print "blast log: "   + output_tmplog

	if run == 1 :
		os.system("makeblastdb  -logfile /dev/null -in "+seq_file+" -dbtype nucl -parse_seqids")
		
		blastcmd='blastn '
		blastcmd=blastcmd+' '
		blastcmd=blastcmd+' -num_descriptions 1000000' 
		blastcmd=blastcmd+' -num_alignments 1000000'
		blastcmd=blastcmd+' -max_target_seqs 1000000'
		blastcmd=blastcmd+' -outfmt 7'
		blastcmd=blastcmd+' -num_threads '   +str(proc)
		blastcmd=blastcmd+' -word_size '   +str(word_size)
		blastcmd=blastcmd+' -query '   +mono_file
		blastcmd=blastcmd+' -db '   +seq_file
		blastcmd=blastcmd+' -out '   +output_tmpblastfile
		blastcmd=blastcmd+' -dust yes '
		blastcmd=blastcmd+' 1>2 > '+output_tmplog
		os.system(blastcmd)
	dt=np.dtype([('query',np.str,100),('target',np.str,100),('simil',np.float), ('length',np.int),('x',np.int), ('y',np.int), ('mono_beg',np.int), ('mono_end',np.int), ('seq_beg', np.int), ('seq_end',np.int), ('pvalue', np.float), ('score', np.float)])
	data = np.genfromtxt(output_tmpblastfile, delimiter='	',
 #    			 names=['monomer','seq','simil', 'length','x', 'y', 'mono_beg', 'mono_end', 'seq_beg', 'seq_end', 'pvalue', 'score']  , 
     			 #dtype=None)
   #  			 dtype=(str,      str,    float, int,      int, int, int,          int,        int,      int,        float,   float))
				 dtype=dt)
#	print data
	if len(np.shape(data)) == 0:
		data=np.asarray([data], dtype=dt)
	if rmout == 1:
		os.remove(output_tmpblastfile)
		os.remove(output_tmplog)
	return data
#############################
def get_seq_length(seq_file):
	if not os.path.exists(seq_file):
    		return -1
	seq_length={}
	for mono_record in SeqIO.parse(seq_file, "fasta"):
		seq_length[mono_record.id]=len(mono_record.seq)
	
	return seq_length

##########################################################
def analyse_sequences_oneByone(work_queue, done_queue, seq_file, mono_file, word_size, limit_lgth_hit, limit_score_hit, limit_distance,verbose=0):

    for seq in iter(work_queue.get, 'STOP'):
	print "#########################################"
	print current_process().name+" : "+seq

	oneseq_file=tempfile.NamedTemporaryFile(delete=False).name
	extract_sequence_to_newFile(seq, seq_file, oneseq_file, name=1)
	
	data=blast(oneseq_file, mono_file, word_size,rmout=1)
	data.sort(order='seq_beg')
	seq_index = SeqIO.index(oneseq_file, "fasta")

	mono_length=get_seq_length(mono_file)
	seq_length=get_seq_length(oneseq_file)

	repeat_region={}
	for hit in data:
		#print hit
	#	print hit['monomer']
		if hit['seq_beg'] < hit['seq_end']:
			str_hit=1
		else:
			tmp=hit['seq_beg']
			hit['seq_beg']=hit['seq_end']
			hit['seq_end']=tmp
			str_hit=-1
		#print str_hit
		lgth_thehit=hit['seq_end']-hit['seq_beg']+1
		
		if lgth_thehit > limit_lgth_hit and hit['score'] > limit_score_hit:
			overlap=-1
			ii=0
 			if hit['seq'] not in repeat_region:
 				repeat_region[hit['seq']]=[]
			while ii < len(repeat_region[hit['seq']]) and overlap == -1:
				if  str_hit == repeat_region[hit['seq']][ii][2] and hit['seq_beg'] <= (repeat_region[hit['seq']][ii][1]+limit_distance) and hit['seq_end'] >= (repeat_region[hit['seq']][ii][0]-limit_distance):
					if hit['seq_beg'] < repeat_region[hit['seq']][ii][0]:
						repeat_region[hit['seq']][ii][0]=hit['seq_beg'] 
					if hit['seq_end'] > repeat_region[hit['seq']][ii][1]:
						repeat_region[hit['seq']][ii][1]=hit['seq_end'] 
					repeat_region[hit['seq']][ii][3]=repeat_region[hit['seq']][ii][1]-repeat_region[hit['seq']][ii][0]+1
					repeat_region[hit['seq']][ii][4]=repeat_region[hit['seq']][ii][4]+1
					overlap=ii
				ii+=1
			
			if overlap == -1:
				repeat_region[hit['seq']].append([hit['seq_beg'], hit['seq_end'],str_hit,hit['seq_end']-hit['seq_beg']+1, 1, hit['seq']])
	
	done_queue.put((seq, repeat_region))
    return True
    
##########################################################  

def writer_to_file(done_queue, seq_file, output_file_seq, output_file_tab, length, verbose=2):
    handle=open(output_file_tab, 'a')
    handle.write("seq begin end strand lgth nbhits\n")
    handle.close()
    seq_index = SeqIO.index(seq_file, "fasta")
    while 1:
        seq, repeat_region = done_queue.get()
	
	if seq != 'STOP':
		handle=open(output_file_tab, 'a')
		handleseq=open(output_file_seq, 'a')
		for chr in repeat_region:
			for rr in repeat_region[chr]:
				if (rr[1]-rr[0]+1) > length:
					handle.write(str(rr[5])+" "+str(rr[0])+" "+str(rr[1])+" "+str(rr[2])+" "+str(rr[3])+" "+str(rr[4])+"\n")

					if rr[2] == -1:
						ss=SeqRecord(seq_index[rr[5]].seq[rr[0]:rr[1]]).reverse_complement().seq
					else:			
						ss=seq_index[rr[5]].seq[rr[0]:rr[1]]
					if args.verbose > 0:
						print(str(rr[5])+" "+str(rr[0])+" "+str(rr[1])+" "+str(rr[2])+" "+str(rr[3])+" "+str(rr[4]))
					
					if args.verbose > 1:
						print(">"+rr[5]+"_"+str(rr[0])+"_"+str(rr[1])+"_"+str(rr[2]))
						print("seq: "+ss)
					handleseq.write(">"+rr[5]+"_"+str(rr[0])+"_"+str(rr[1])+"_"+str(rr[2])+"\n")
					handleseq.write(str(ss)+"\n")
			
		handle.close()
		handleseq.close()
	else:
		break
    return True


#############################
parser = argparse.ArgumentParser(description='Search for tandemly repeated regions ...', 
epilog='This program searchs for repeated regions including occurences of monomers. Searched sequences (-s option) are compared with monomers (-m option) by using blast (needed).')


parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-m', dest="mono_file", required=True)
parser.add_argument('-t', dest="output_tabfile", required=True)
parser.add_argument('-o', dest="output_seqfile", required=True)
parser.add_argument('-p',  dest="proc",   type=int, default=4)

parser.add_argument('-W',  dest="word_size",   type=int, default=4)
parser.add_argument('-d',  dest="distance",    type=int, default=100)
parser.add_argument('-S',  dest="score_hit",   type=float, default=0.0)
parser.add_argument('-L',  dest="lgth_hit",   type=int, default=0)
parser.add_argument('-l',  dest="lgth_seq",   type=int, default=100)

parser.add_argument('-v',  dest="verbose",   type=int, default=0)

args = parser.parse_args()

if os.path.exists(args.output_tabfile):
    os.remove(args.output_tabfile)
if os.path.exists(args.output_seqfile):
    os.remove(args.output_seqfile)

data=blast(args.seq_file, args.mono_file, proc=args.proc, word_size=args.word_size, output_tmpblastfile="", rmout=0)
print data

#data.sort(order='target')
#seq_index = SeqIO.index(oneseq_file, "fasta")

mono_length=get_seq_length(args.mono_file)
seq_length=get_seq_length(args.seq_file)

print seq_length
repeat_region={}
for hit in data:
	print hit
	if hit['seq_beg'] < hit['seq_end']:
		str_hit=1
	else:
		tmp=hit['seq_beg']
		hit['seq_beg']=hit['seq_end']
		hit['seq_end']=tmp
		str_hit=-1
	lgth_thehit=hit['seq_end']-hit['seq_beg']+1
		
	if lgth_thehit > args.lgth_hit and hit['score'] > args.score_hit:
		overlap=-1
		ii=0
 		if hit['target'] not in repeat_region:
 			repeat_region[hit['target']]=[]
		while ii < len(repeat_region[hit['target']]) and overlap == -1:
			if  str_hit == repeat_region[hit['target']][ii][2] and hit['seq_beg'] <= (repeat_region[hit['target']][ii][1]+args.distance) and hit['seq_end'] >= (repeat_region[hit['target']][ii][0]-args.distance):
				if hit['seq_beg'] < repeat_region[hit['target']][ii][0]:
					repeat_region[hit['target']][ii][0]=hit['seq_beg'] 
				if hit['seq_end'] > repeat_region[hit['target']][ii][1]:
					repeat_region[hit['target']][ii][1]=hit['seq_end'] 
				repeat_region[hit['target']][ii][3]=repeat_region[hit['target']][ii][1]-repeat_region[hit['target']][ii][0]+1
				repeat_region[hit['target']][ii][4]=repeat_region[hit['target']][ii][4]+1
				overlap=ii
			ii+=1
			
		if overlap == -1:
			repeat_region[hit['target']].append([hit['seq_beg'], hit['seq_end'],str_hit,hit['seq_end']-hit['seq_beg']+1, 1, hit['query']])



handle=open(args.output_tabfile, 'a')
handleseq=open(args.output_seqfile, 'a')
handle.write("seq begin end strand lgth nbhits\n")
seq_index = SeqIO.index(args.seq_file, "fasta")
for chr in repeat_region:
			for rr in repeat_region[chr]:
				if (rr[1]-rr[0]+1) > args.lgth_seq:
					handle.write(str(chr)+" "+str(rr[0])+" "+str(rr[1])+" "+str(rr[2])+" "+str(rr[3])+" "+str(rr[4])+"\n")

					if rr[2] == -1:
						
						ss=SeqRecord(seq_index[chr].seq[rr[0]:rr[1]]).reverse_complement().seq
					else:			
						ss=seq_index[chr].seq[rr[0]:rr[1]]
					if args.verbose > 0:
						print(str(chr)+" "+str(rr[0])+" "+str(rr[1])+" "+str(rr[2])+" "+str(rr[3])+" "+str(rr[4]))
					
					if args.verbose > 1:
						print(">"+chr+"_"+str(rr[0])+"_"+str(rr[1])+"_"+str(rr[2]))
						print("seq: "+ss)
					handleseq.write(">"+chr+"_"+str(rr[0])+"_"+str(rr[1])+"_"+str(rr[2])+"\n")
					handleseq.write(str(ss)+"\n")
			
handle.close()
handleseq.close()





print "Bye bye  :-) !!!"

