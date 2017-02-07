#! /usr/bin/env python

import re
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
iport my_function

def blast(seq_file, mono_file, proc, word_size, output_tmpblastfile="", rmout=1):
	#############################
	# blasting ...
	print 'blasting ...'
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
		os.system("formatdb -i "+seq_file+" -p F -o")
		blastcmd='blastall -p blastn '
		blastcmd=blastcmd+' -F F'
		blastcmd=blastcmd+' -v 100000' 
		blastcmd=blastcmd+' -b 100000'
		blastcmd=blastcmd+' -m 8'
		blastcmd=blastcmd+' -a '   +str(proc)
		blastcmd=blastcmd+' -W '   +str(word_size)
		blastcmd=blastcmd+' -i '   +mono_file
		blastcmd=blastcmd+' -d '   +seq_file
		blastcmd=blastcmd+' -o '   +output_tmpblastfile
		blastcmd=blastcmd+' 1>2 > '+output_tmplog
		os.system(blastcmd)
	
	data = np.genfromtxt(output_tmpblastfile, delimiter='	',
     			 names=['monomer','seq','simil', 'length','x', 'y', 'mono_beg', 'mono_end', 'seq_beg', 'seq_end', 'pvalue', 'score']  , dtype=None)
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

#############################
parser = argparse.ArgumentParser(description='Search for tandemly repeated regions ...')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-m', dest="mono_file", required=True)
parser.add_argument('-t', dest="output_tabfile", required=True)
parser.add_argument('-o', dest="output_seqfile", required=True)
parser.add_argument('-B', dest="blastfile", default="")
parser.add_argument('-r',  dest="rm_blastfile",   type=int, default=1)
parser.add_argument('-W',  dest="word_size",   type=int, default=4)
parser.add_argument('-p',  dest="processors",  type=int, default=4)
parser.add_argument('-d',  dest="distance",    type=int, default=100)
parser.add_argument('-S',  dest="score_hit",   type=float, default=0.0)
parser.add_argument('-L',  dest="lgth_hit",   type=int, default=0)
parser.add_argument('-v',  dest="verbose",   type=int, default=0)

args = parser.parse_args()
output_tmpblastfile=tempfile.NamedTemporaryFile(delete=False).name




if os.path.exists(args.output_tabfile):
    os.remove(args.output_tabfile)
if os.path.exists(args.output_seqfile):
    os.remove(args.output_seqfile)

data=blast(args.seq_file, args.mono_file, args.processors, args.word_size,args.blastfile, args.rm_blastfile)
if data == -1:
	print 'Sequence file '+args.seq_file+' doesn\'t exist (option -s) !'
	quit()
if data == -2:
	print 'Monomer file '+args.mono_file+' doesn\'t exist (option -m) !'
	quit()
	
data.sort(order='seq_beg')
seq_index = SeqIO.index(args.seq_file, "fasta")

mono_length=get_seq_length(args.mono_file)
seq_length=get_seq_length(args.seq_file)

repeat_region={}
for hit in data:


#	for ss in repeat_region:
#		print ss, len(repeat_region[ss])

	if hit['seq_beg'] < hit['seq_end']:
		str_hit=1
	else:
		tmp=hit['seq_beg']
		hit['seq_beg']=hit['seq_end']
		hit['seq_end']=tmp
		str_hit=-1

	lgth_hit=hit['seq_end']-hit['seq_beg']+1
		
	if lgth_hit > args.lgth_hit and hit['score'] > args.score_hit:
		overlap=-1
		ii=0
 		if hit['seq'] not in repeat_region:
 			repeat_region[hit['seq']]=[]
		while ii < len(repeat_region[hit['seq']]) and overlap == -1:
			if  str_hit == repeat_region[hit['seq']][ii][2] and hit['seq_beg'] <= (repeat_region[hit['seq']][ii][1]+args.distance) and hit['seq_end'] >= (repeat_region[hit['seq']][ii][0]-args.distance):
				if hit['seq_beg'] < repeat_region[hit['seq']][ii][0]:
					repeat_region[hit['seq']][ii][0]=hit['seq_beg'] 
				if hit['seq_end'] > repeat_region[hit['seq']][ii][1]:
					repeat_region[hit['seq']][ii][1]=hit['seq_end'] 
				repeat_region[hit['seq']][ii][3]=repeat_region[hit['seq']][ii][1]-repeat_region[hit['seq']][ii][0]+1
				repeat_region[hit['seq']][ii][4]=repeat_region[hit['seq']][ii][4]+1
 #				print repeat_region[hit['seq']][ii]
				overlap=ii
			ii+=1
			
		if overlap == -1:
			repeat_region[hit['seq']].append([hit['seq_beg'], hit['seq_end'],str_hit,hit['seq_end']-hit['seq_beg']+1, 1, hit['seq']])
#	print hit
handle=open(args.output_tabfile,'w')
handleseq=open(args.output_seqfile,'w')

handle.write("seq begin end strand lgth nbhits\n")

for chr in repeat_region:
	
	print "###########"
	print chr
	
	mmL=mmN=MML=MMN=-1
	if len(repeat_region[chr]) >= 1:
		mmL=repeat_region[chr][0][3]
		MML=repeat_region[chr][0][3]
		mmN=repeat_region[chr][0][4]
		MMN=repeat_region[chr][0][4]
		
	for rr in repeat_region[chr]:
		if rr[3] < mmL:
			mmL=rr[3]
		if rr[3] > MML:
			MML=rr[3]
		if rr[4] < mmN:
			mmN=rr[4]
		if rr[4] > MMN:
			MMN=rr[4]
			
		handle.write(str(rr[5])+" "+str(rr[0])+" "+str(rr[1])+" "+str(rr[2])+" "+str(rr[3])+" "+str(rr[4])+"\n")
	
	
	
	#print "##########################"
	#print chr, 
	#print seq_length[chr], len(repeat_region[chr]), mmL, MML, mmN, MMN
	
	for rr in repeat_region[chr]:
		if rr[2] == -1:
			ss=SeqRecord(seq_index[rr[5]].seq[rr[0]:rr[1]]).reverse_complement().seq
		else:			
			ss=seq_index[rr[5]].seq[rr[0]:rr[1]]
		if args.verbose > 0:
			print(">"+rr[5]+"_"+str(rr[0])+"_"+str(rr[1])+"_"+str(rr[2]))
			print("seq"+ss)
		handleseq.write(">"+rr[5]+"_"+str(rr[0])+"_"+str(rr[1])+"_"+str(rr[2])+"\n")
		handleseq.write(str(ss)+"\n")
		
handle.close()
handleseq.close()
print "Bye bye :-) !!"
