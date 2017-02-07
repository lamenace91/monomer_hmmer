#! /usr/bin/env python

from scipy import stats
import csv
import os
import sys
import re
import numpy as np
from numpy import genfromtxt
from Bio import AlignIO
import time
import argparse
import tempfile
from Bio import SeqIO
import math
#from scipy import DataFrame
from multiprocessing import Lock, Process, Queue, current_process

#############################
parser = argparse.ArgumentParser(description='pairwise alignment of monomers ...')
parser.add_argument('-d', dest="dist_file")
parser.add_argument('-s', dest="seq_file")
args = parser.parse_args()

sequence_dict = SeqIO.index(args.seq_file, "fasta")

bloc_dist={}
with open(args.dist_file,'r') as f:
#with open("zzz",'r') as f:
    reader=csv.reader(f, delimiter=' ')
    nb=0
    ss=''
    for row in reader:
    	bloc=row[0].partition(":")[0].replace('lcl|','')
	
    	if bloc not in bloc_dist:
		bloc_dist[bloc]=list()		
	bloc_dist[bloc].append(float(row[4]))


bloc_dict={}
for ss in sequence_dict.keys():
	bloc=ss.partition(":")[0].replace('lcl|','')
	if bloc not in bloc_dict:
		bloc_dict[bloc]=[bloc, 0, -1.0, re.split("[_]", ss)[0].replace('lcl|',''),    int(re.split("[_]", bloc)[-3]),     int(re.split("[_]", bloc)[-2]),   int(re.split("[_]", ss)[-2])-int(re.split("[_]", bloc)[-3])    ]
		
		if bloc in bloc_dist:
#			print "###########################"
#			print bloc
#			print bloc_dict[bloc]
#			print bloc_dict[bloc][2]
#			print len(bloc_dist[bloc])
			#print bloc_dist[bloc]
#			print np.mean(bloc_dist[bloc])
			bloc_dict[bloc][2]=np.mean(bloc_dist[bloc])
	bloc_dict[bloc][1]=bloc_dict[bloc][1]+1
	

for bb in bloc_dict:
	print bloc_dict[bb]


print "Bye bye !!!"
quit()
