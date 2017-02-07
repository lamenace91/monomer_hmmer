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
 
def gen_sub(s, len_chunk):
   	for start in range(0, len(s)-len_chunk+1):
       		yield s[start:start+len_chunk]
   #############################
def get_distance(wc_seq1, wc_seq2):
	dd=0
	for kmer in wc_seq1:
		dd=dd+(wc_seq1[kmer]-wc_seq2[kmer])*(wc_seq1[kmer]-wc_seq2[kmer])
	return(dd)

#############################
#############################
parser = argparse.ArgumentParser(description='K-mers count calculation.')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-o', dest="output_file", required=True)
parser.add_argument('-k', dest="k",   type=int, default=2)
parser.add_argument('-f', dest="freq",   action='store_true', default=False)
parser.add_argument('-v', dest="verbose",   type=int, default=1)



N=20

args = parser.parse_args()

if os.path.exists(args.output_file):
    os.remove(args.output_file)

handle1=open(args.seq_file, 'r')
seq_dict = SeqIO.to_dict(SeqIO.parse(handle1, "fasta"))
handle1.close()

if args.verbose > 0:
	print "sequences: ", len(seq_dict)




#kmers={}
#for ii in range(args.k1,args.k2+1): 
#	kmers[ii]=list(product('ACGT', repeat=ii))
#	print str(ii)+"-mers: ", len(kmers[ii])
#	
#	for jj in range(0, len(kmers[ii])):
#		kmers[ii][jj]=''.join(kmers[ii][jj])
	
all_words=[]
kmers_count={}
for ss1 in seq_dict:
	kmers_count[ss1]=Counter([sub for sub in gen_sub(str(seq_dict[ss1].seq.upper()),args.k )])
	all_words= all_words+ kmers_count[ss1].keys()
	all_words=list(set(all_words))
	if args.freq:
		for ii in kmers_count[ss1]:
			kmers_count[ss1][ii]=float(kmers_count[ss1][ii])/(len(seq_dict[ss1].seq)-args.k+1.0)
#		print ss1, ii, kmers_count[ss1][ii]
		
#		print "#######################"
#		print ss1
#		print kmers_count[ss1]
#		print "#######################"
#	print all_words
		
#ss1='lcl|chr12_38292356_38301939_C:7635-7814'
#ss2='lcl|chr12_38292356_38301939_C:6265-6435'
#print kmers_count[ss1]-kmers_count[ss2]	
#xx=kmers_count[ss1]
#xx.subtract(kmers_count[ss2])
#print xx
#print "#######################"
all_words=sorted(all_words)

handle=open(args.output_file,'w')
handle.write("id")
for kkk in all_words:
	handle.write (" "+kkk)
handle.write ("\n")
	
for ss1 in seq_dict:
	handle.write(ss1)
	for kkk in all_words:
		if kkk in kmers_count[ss1]:
			handle.write (" "+str(kmers_count[ss1][kkk]))
		else:
			if args.freq:
				handle.write (" 0.0")
			else:
				handle.write (" 0")
	handle.write ("\n")

handle.close()
print "\nBye bye !!!"

