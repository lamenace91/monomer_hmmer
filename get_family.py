#! /usr/bin/env python

from collections import Counter
import csv
import os
import sys
import re
import numpy as np
from itertools import *
from sys import stdout
from numpy import genfromtxt
from Bio import AlignIO
import time
import argparse
import tempfile
from Bio import SeqIO
import math
#from scipy import DataFrame
from multiprocessing import Lock, Process, Queue, current_process
from progress_bar import ProgressBar
#############################

#############################################################
def get_similarity_and_score_length(seq1, seq2):
	align_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq1_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq2_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	count = SeqIO.write(seq1, seq1_tmpfile, "fasta")
	count = SeqIO.write(seq2, seq2_tmpfile, "fasta")
	os.system("needle  -asequence "+seq1_tmpfile+" -bsequence "+seq2_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 0.5 -aformat3 fasta 2>1 > /dev/null")
	f=open(align_tmpfile,"r")
	alignment = AlignIO.read(align_tmpfile, "fasta")
	f.close()
	simil=compute_similarity(alignment[0].seq,alignment[1].seq)
	score=compute_score(alignment[0].seq,alignment[1].seq)

	os.remove(align_tmpfile)
	os.remove(seq1_tmpfile)
	os.remove(seq2_tmpfile)
	return (simil,score, len(alignment[0].seq))

#############################################################
def compute_score(seq1, seq2):
	go=-1
	ge=-0.5
	score=0
	char={'A':0,'T':1,'G':2,'C':3,'S':4,'W':5,'R':6,'Y':7,'K':8,'M':9,'B':10,'V':11,'H':12,'D':13,'N':14}
 	dist=[
[  5,  -4,  -4,  -4,  -4,   1,   1,  -4,  -4,   1,  -4,  -1,  -1,  -1,  -2],
[ -4,   5,  -4,  -4,  -4,   1,  -4,   1,   1,  -4,  -1,  -4,  -1,  -1,  -2],
[ -4,  -4,   5,  -4,   1,  -4,   1,  -4,   1,  -4,  -1,  -1,  -4,  -1,  -2],
[ -4,  -4,  -4,   5,   1,  -4,  -4,   1,  -4,   1,  -1,  -1,  -1,  -4,  -2],
[ -4,  -4,   1,   1,  -1,  -4,  -2,  -2,  -2,  -2,  -1,  -1,  -3,  -3,  -1],
[  1,   1,  -4,  -4,  -4,  -1,  -2,  -2,  -2,  -2,  -3,  -3,  -1,  -1,  -1],
[  1,  -4,   1,  -4,  -2,  -2,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -1,  -1],
[ -4,   1,  -4,   1,  -2,  -2,  -4,  -1,  -2,  -2,  -1,  -3,  -1,  -3,  -1],
[ -4,   1,   1,  -4,  -2,  -2,  -2,  -2,  -1,  -4,  -1,  -3,  -3,  -1,  -1],
[  1,  -4,  -4,   1,  -2,  -2,  -2,  -2,  -4,  -1,  -3,  -1,  -1,  -3,  -1],
[ -4,  -1,  -1,  -1,  -1,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -2,  -2,  -1],
[ -1,  -4,  -1,  -1,  -1,  -3,  -1,  -3,  -3,  -1,  -2,  -1,  -2,  -2,  -1],
[ -1,  -1,  -4,  -1,  -3,  -1,  -3,  -1,  -3,  -1,  -2,  -2,  -1,  -2,  -1], 
[ -1,  -1,  -1,  -4,  -3,  -1,  -1,  -3,  -1,  -3,  -2,  -2,  -2,  -1,  -1],
[ -2,  -2,  -2,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1]]
	
	ii=0
	seq1=seq1.upper()
	seq2=seq2.upper()	
	while ii < len(seq1):
		if seq1[ii] == '-' and ii==0:
			score+=go
		elif seq2[ii] == '-' and ii==0:
			score+=go
		elif seq1[ii] == '-' and seq1[ii-1] != '-':
			score+=go
		elif seq2[ii] == '-' and seq2[ii-1] != '-':
			score+=go
		elif seq1[ii] == '-' and seq1[ii-1] == '-':
			score+=ge
		elif seq2[ii] == '-' and seq2[ii-1] == '-':
			score+=ge
		else:
			score+=dist[char[seq1[ii]]][char[seq2[ii]]]
#		print str(ii), seq1[ii], seq2[ii], str(score)
		ii+=1
	return score
#############################################################
def compute_similarity(seq1, seq2):
	simil=0.0
	total=0.0
	ii=0
	seq1=seq1.upper()
	seq2=seq2.upper()	
	while ii < len(seq1):
		if seq1[ii] != '-' and seq2[ii] != '-':
			if seq1[ii] == seq2[ii]:
				simil+=1.0
			total+=1.0
		ii+=1
	if total > 0:
		return simil/total
	else:
		return -1.0


##########################################
def writer_to_file(done_queue, output_file):    

    while 1:
        m = done_queue.get()
  	handle1=open(output_file, 'a')
	handle1.write(str(m)+'\n')
   	handle1.close()
    return True

#############################
def worker1(work_queue, done_queue,mono_dict):
    for ss in iter(work_queue.get, 'STOP'):
#	print "#########################################"
#	print current_process().name+" : "+ss.id
	result=[ss.id]
	for mono in mono_dict:
		(simil,score, length)=get_similarity_and_score_length(ss,mono_dict[mono])
		print ss.id+" "+mono+" "+str(simil)
		result.append((mono, simil, score))
	done_queue.put(result)
    return True




#############################
parser = argparse.ArgumentParser(description='pairwise alignment of monomers ...')
parser.add_argument('-s', dest="seq_file")
parser.add_argument('-m', dest="mono_file")
parser.add_argument('-o', dest="output_file")
parser.add_argument('-p', dest="proc",   type=int, default=8)

args = parser.parse_args()

if os.path.exists(args.output_file):
    os.remove(args.output_file)

handle1=open(args.seq_file, 'r')
seq_dict = SeqIO.to_dict(SeqIO.parse(handle1, "fasta"))
handle1.close()


handle1=open(args.mono_file, 'r')
mono_dict = SeqIO.to_dict(SeqIO.parse(handle1, "fasta"))
handle1.close()

workers = args.proc
work_queue = Queue()
done_queue = Queue()
processes = []


q = Process(target=writer_to_file, args=(done_queue,args.output_file))
q.start()

for ss in seq_dict.keys():
	work_queue.put(seq_dict[ss])	


print "seq: "+str(len(seq_dict))
print "mono: "+str(len(mono_dict))

for w in xrange(workers):
        p = Process(target=worker1, args=(work_queue, done_queue, mono_dict) )
        p.start()
        processes.append(p)
        work_queue.put('STOP')
	
for p in processes:
        p.join()

done_queue.put('STOP')
		
print "Bye bye !!!"


