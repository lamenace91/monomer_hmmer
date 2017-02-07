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
parser = argparse.ArgumentParser(description='K-mers distance calculation.')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-o', dest="output_file", required=True)
parser.add_argument('-k', dest="k",   type=int, default=2)
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

bloc_dict={}

for ss in seq_dict.keys():
	bloc=ss.partition(":")[0]
	if bloc not in bloc_dict:
		bloc_dict[bloc]=list()
	bloc_dict[bloc].append(ss)
	
for bb in bloc_dict:
	bloc_dict[bb].sort(key= lambda tup:tup[1])


#kmers={}
#for ii in range(args.k1,args.k2+1): 
#	kmers[ii]=list(product('ACGT', repeat=ii))
#	print str(ii)+"-mers: ", len(kmers[ii])
#	
#	for jj in range(0, len(kmers[ii])):
#		kmers[ii][jj]=''.join(kmers[ii][jj])
	

kmers_count={}
for bb in bloc_dict:
	for ss1 in bloc_dict[bb]:
		kmers_count[ss1]=Counter([sub for sub in gen_sub(str(seq_dict[ss1].seq),args.k )])
		for ii in kmers_count[ss1]:
			kmers_count[ss1][ii]=float(kmers_count[ss1][ii])/(len(seq_dict[ss1].seq)-args.k+1.0)
#		print "#######################"
#		print ss1
#		print kmers_count[ss1]
#		print "#######################"
		
#ss1='lcl|chr12_38292356_38301939_C:7635-7814'
#ss2='lcl|chr12_38292356_38301939_C:6265-6435'
#print kmers_count[ss1]
#print kmers_count[ss2]
#print kmers_count[ss1]-kmers_count[ss2]	
#xx=kmers_count[ss1]
#xx.subtract(kmers_count[ss2])
#print xx
#print "#######################"
nb=0.0
handle=open(args.output_file,'w')
nb_dd=len(seq_dict)*len(seq_dict)
if args.verbose > 0:
	print "nb of pairs: ", nb_dd
handle.write("id seq1 seq2"+" "+str(args.k)+"-mers\n")
p = ProgressBar(nb_dd)	
progresspc=-1
for bb1 in bloc_dict:
	for bb2 in bloc_dict:
		for ss1 in bloc_dict[bb1]:
			for ss2 in bloc_dict[bb2]:
				handle.write(str(nb)+" "+ss1+" "+ss2)
				nb+=1.0
				if args.verbose > 0:
					#print int(nb/nb_dd*100.0), progresspc
					if int(nb/nb_dd*100.0) > progresspc:
						progresspc=int(nb/nb_dd*100.0)
						progress=int(nb/nb_dd*N)
						print '\r[{0}{1}] {2}%'.format('#'*(progress),' '*(N-progress) , progresspc), 
						sys.stdout.flush()	
				kmers_min=kmers_count[ss2]|kmers_count[ss1]
				kmers_max=kmers_count[ss2]&kmers_count[ss1]
				
				dd=0
				for ii in kmers_min:
					dd=dd+(kmers_max[ii]-kmers_min[ii])*(kmers_max[ii]-kmers_min[ii])
				dd1=0
				sum1=kmers_count[ss2]-kmers_count[ss1]
				sum2=kmers_count[ss1]-kmers_count[ss2]
				for zz in sum1:
					dd1=dd1+sum1[zz]*sum1[zz]
				for zz in sum2:
					dd1=dd1+sum2[zz]*sum2[zz]
				if args.verbose > 1:
					print ss1, ss2, dd, dd1
				handle.write(" "+str(dd))
				handle.write("\n")
handle.close()
print "\nBye bye !!!"
