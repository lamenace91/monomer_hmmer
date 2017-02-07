#! /usr/bin/env python

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
    #############################
def get_distance(wc_seq1, wc_seq2):
	dd=0
	for kmer in wc_seq1:
		dd=dd+(wc_seq1[kmer]-wc_seq2[kmer])*(wc_seq1[kmer]-wc_seq2[kmer])
	return(dd)

#############################
#############################
parser = argparse.ArgumentParser(description='pairwise alignment of monomers ...')
parser.add_argument('-s', dest="seq_file")
parser.add_argument('-o', dest="output_file")
parser.add_argument('-k1', dest="k1",   type=int, default=2)
parser.add_argument('-k2', dest="k2",   type=int, default=5)


args = parser.parse_args()

if os.path.exists(args.output_file):
    os.remove(args.output_file)

handle1=open(args.seq_file, 'r')
seq_dict = SeqIO.to_dict(SeqIO.parse(handle1, "fasta"))
handle1.close()


print "sequences: ", len(seq_dict)

bloc_dict={}

for ss in seq_dict.keys():
	bloc=ss.partition(":")[0]
	if bloc not in bloc_dict:
		bloc_dict[bloc]=list()
	bloc_dict[bloc].append(ss)
	
for bb in bloc_dict:
	bloc_dict[bb].sort(key= lambda tup:tup[1])


kmers={}
for ii in range(args.k1,args.k2+1): 
	kmers[ii]=list(product('ACGT', repeat=ii))
	print str(ii)+"-mers: ", len(kmers[ii])
	
	for jj in range(0, len(kmers[ii])):
		kmers[ii][jj]=''.join(kmers[ii][jj])
	

kmers_count={}
for bb in bloc_dict:
	for ss1 in bloc_dict[bb]:
		kmers_count[ss1]={}
		for ii in range(args.k1,args.k2+1): 
			kmers_count[ss1][ii]={}
			for kmer in range (0,len(kmers[ii])):
				#print kmer, seq_dict[ss1].seq.count(kmer)
				count function don't take account overlapping strings
				kmers_count[ss1][ii][kmers[ii][kmer]]=float(seq_dict[ss1].seq.count(kmers[ii][kmer]))/(len(seq_dict[ss1].seq)-ii+1.0)
				print ss1, kmers[ii][kmer], kmers_count[ss1][ii][kmers[ii][kmer]], seq_dict[ss1].seq.count(kmers[ii][kmer])
#				print kmers_count[ss1][ii][kmer]
				
nb=0
handle=open(args.output_file,'w')
nb_dd=len(seq_dict)*len(seq_dict)
handle.write("seq begin end strand lgth nbhits")
	handle.write(" "+str(ii)+"-mers")
handle.write("\n")



p = ProgressBar(nb_dd)	
for bb1 in bloc_dict:
#	for bb2 in bloc_dict:
		for ss1 in bloc_dict[bb1]:
#			for ss2 in bloc_dict[bb2]:
			for ss2 in bloc_dict[bb1]:
				handle.write(str(nb)+" "+ss1+" "+ss2)
				for ii in range(args.k1,args.k2+1): 
					nb+=1
					if (nb%(nb_dd/1000000.0))==0:
						p.update_time(nb)
						stdout.write("\r%s" % p)
					#update_progress(30*nb/nb_dd)
					print "###############"
					print kmers_count[ss1][ii]
					print kmers_count[ss2][ii]
					
					#print [(kmers_count[ss1][ii][j]-kmers_count[ss2][ii][i])*(kmers_count[ss1][ii][j]-kmers_count[ss2][ii][i]) for i,j in zip(kmers_count[ss1][ii],kmers_count[ss2][ii])]
					dd=sum([(kmers_count[ss1][ii][j]-kmers_count[ss2][ii][i])*(kmers_count[ss1][ii][j]-kmers_count[ss2][ii][i]) for i,j in zip(kmers_count[ss1][ii],kmers_count[ss2][ii])])
					print ss1, ss2, dd
					handle.write(" "+str(dd))
				handle.write("\n")
handle.close()
print "Bye bye !!!"
