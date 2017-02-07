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
import operator
####################################################################
def ratio_distances(hor_l, hor_beg, distances):
	
	dist_intra=[]
	dist_inter=[]
	selected=[]
	for ii in range(0,hor_l):
		for jj in hor_beg:
			for kk in hor_beg:
				if jj != kk and jj < kk:
					selected.append((kk+ii,jj+ii))
					selected.append((jj+ii,kk+ii))
					if (jj+ii) < len(distances[0]) and (kk+ii) < len(distances[0]):	
#						print len(distances[0]), ii, jj+ii, kk+ii
						print kk, jj, ii, distances[jj+ii][kk+ii]
						dist_intra.append(distances[jj+ii][kk+ii])
	for jj in range(0,len(distances[0])):
		for kk in range(0,len(distances[0])):
			if jj != kk and jj < kk and (jj,kk) not in selected :
				dist_inter.append(distances[jj][kk])
		
#	print "intra", dist_intra
#	print "############################"
#	print "inter", dist_inter
	tt=(0,1)
	if len(dist_intra) > 1 and len(dist_inter) > 1:
		tt=stats.ttest_ind(dist_intra, dist_inter)
	
	return 	len(dist_intra), np.mean(dist_intra), len(dist_inter),np.mean(dist_inter),np.mean(dist_intra)/np.mean(dist_inter), tt[0], tt[1]

####################################################################



parser = argparse.ArgumentParser(description='search HOR from a distance file. Distance have to be include at least five columns, with seqnames in first ones and similarities in fifth. See for example the output of align2monomers.py.')
#parser.add_argument('-s', dest="seq_file", required=true)
parser.add_argument('-d', dest="dist_file", required=True)
#parser.add_argument('-i', dest="info_file")
#parser.add_argument('-a', dest="all_blocks", default=True)
parser.add_argument('-o', dest="output_file", required=True)
args = parser.parse_args()


bloc_dict={}
posseq_dict={}
with open(args.dist_file,'r') as f:
#with open("zzz",'r') as f:
    reader=csv.reader(f, delimiter=' ')
    nb=0
    ss=''
    for row in reader:
    	bloc1=row[0].partition(":")[0].replace('lcl|','')
    	bloc2=row[1].partition(":")[0].replace('lcl|','')
	if bloc1 == bloc2:
		bloc=bloc1
		seq1=row[0].replace('lcl|','')
		seq2=row[1].replace('lcl|','')
		print seq1, seq2, bloc, row[0][1], row[1][1]
    		if bloc not in bloc_dict:
			bloc_dict[bloc]=list()
			posseq_dict[bloc]=dict()
  		
		bloc_dict[bloc].append((seq1, seq2, float(row[4]), int(re.split("[-:]", row[0])[1]),int(re.split("[-:]", row[1])[1]) ))
		posseq_dict[bloc][seq1]=int(re.split("[-:]", row[0])[1])
		posseq_dict[bloc][seq2]=int(re.split("[-:]", row[1])[1])
#		print row[0].replace('lcl|','')
		#print bloc_dict[bloc]
f.close()
print len(bloc_dict)
handle=open(args.output_file, 'w')

orderseq_dict={}
matrix_dict={}

for bb in bloc_dict:
	orderseq_dict[bb]=dict()
	matrix_dict[bb]=dict()

	ll=len(posseq_dict[bb])
#	print ll
#	print posseq_dict[bb]
#	print "#################"
	sorted_seq = sorted(posseq_dict[bb].iteritems(), key=operator.itemgetter(1))
	print bb, posseq_dict[bb]

	print bb, sorted_seq
	for ii in range(0,ll):
#		print ii, "---->", sorted_seq[ii][0]
		orderseq_dict[bb][sorted_seq[ii][0]]=ii
#	print bb, sorted_seq
#	print bb, orderseq_dict[bb]

	matrix_dict[bb] = np.zeros((ll,ll))
	for dd in bloc_dict[bb]:
#		print dd[0], orderseq_dict[bb][dd[0]]
#		print dd[1], orderseq_dict[bb][dd[1]]
#		print dd[2]
		matrix_dict[bb][orderseq_dict[bb][dd[0]]][orderseq_dict[bb][dd[1]]]=dd[2]
#	print matrix_dict[bb]
	head="array nb_mono hor_order q simil_hor nb_simil simil tvalue pvalue\n" 
	handle.write(head)
	for ii in range(2,min(40, ll)):
		print ii
		out=bb+" "+str(ll)+" "+str(ii)
		rd=ratio_distances(ii, range(0,ll,ii),matrix_dict[bb] )
		for ii in rd:
			out=out+" "+str(ii)
		handle.write(out+'\n')
		handle.flush()
		print out
handle.close()
