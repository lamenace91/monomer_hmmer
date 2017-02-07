#! /usr/bin/env python

import sys
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

#############################################################
## selected non-overlapping hits
#############################################################
def get_nonoverlapping_monomers(data, blast_maximal_overlap):

	nb_selected_monomers = 0
	selected_monomers=[]
	if len(data) == 0:
		return data
	ind = np.argmax(data[0]['score'])
	while len(selected_monomers) < len(data) and data[ind]['score'] > 0:
		overlap=0
		for sm in  selected_monomers:
			if data[ind]['seq_beg'] < (data[sm]['seq_end']-blast_maximal_overlap) and data[ind]['seq_end'] >  (data[sm]['seq_beg']+blast_maximal_overlap):
				overlap=1
		if overlap == 0:
			nb_selected_monomers += 1
			selected_monomers.append(ind)
			if nb_selected_monomers == 1:
				selected_monomers_data=data[ind]
				selected_monomers_data=np.append(selected_monomers_data,data[ind])
				selected_monomers_data=np.delete(selected_monomers_data, 1, 0)
			else:
				selected_monomers_data= np.append(selected_monomers_data,data[ind])
		data[ind]['score']=-1
		if len(data) == 1:
			ind=0
		else:
			ind = np.argmax(data[:]['score'])	
	if len(selected_monomers_data) > 1:
		selected_monomers_data.sort(order='seq_beg')
		
	return(selected_monomers_data)

######################################
## Format a BLAST db 
######################################
def formatdb_nuc(seq_file):
	os.system("formatdb -i "+seq_file+" -p F -o")
	
######################################
## Run a blast and parse the output
######################################
def blastN_and_parseOutput(seq_file, database, options="-F F -v 100000 -b 100000 -m 8"):
	blast_file =tempfile.NamedTemporaryFile(delete=False).name
	log_file    =tempfile.NamedTemporaryFile(delete=False).name
	
	blastcmd='blastall -p blastn '
	blastcmd=blastcmd          + options
	blastcmd=blastcmd+' -i '   +seq_file
	blastcmd=blastcmd+' -d '   +database
	blastcmd=blastcmd+' -o '   +blast_file
	blastcmd=blastcmd+' 1>2 > '+ log_file
	os.system(blastcmd)
	
      	dt=[('seq',object) ,('monomer',object), ('simil','f'), ('length','i') ,('mono_beg','i'), ('mono_end','i'), ('seq_beg','i'), ('seq_end','i'), ('evalue','f'),('score','f')]
	data=np.array([], dtype=dt)
	if os.stat(blast_file)[6] > 0:
		data = np.genfromtxt(blast_file, delimiter='\t', usecols=[0,1,2,3,6,7,8,9,10,11],
      			dtype=dt)
	os.remove(blast_file)
	os.remove(log_file)
	return data

######################################
## Parse monomers/blast_hits for minimal values
######################################
def select_monomers_on_length(data, length):
	data=data[data['length'] >= length]
	return data
	
def select_monomers_on_score(data, score):
	data=data[data['score'] >= score]
	return data
	
def select_monomers_on_evalue(data, evalue):
	data=data[data['evalue'] <= evalue]
	return data


######################################
## Extract a sequence into a new file
######################################
def extract_sequence_to_newFile(seq, from_file, to_file):
	found=0
	record_index = SeqIO.index(from_file, "fasta")
	if seq in record_index:
		output_handle = open(to_file, "w")
		SeqIO.write(record_index[seq] , output_handle, "fasta")
		output_handle.close()		
		return 1
	else:
		return 0	
	return found

######################################
## Run analysis for one sequence
######################################
def sequence_analysis(seq_name, seq_file, mono_file, blast_maximal_overlap, blast_minimal_length, blast_minimal_score, blast_maximal_evalue, waterman_minimal_score, waterman_minimal_length):
	print seq_name
	
	this_seq_file =tempfile.NamedTemporaryFile(delete=False).name
	if extract_sequence_to_newFile(seq_name,seq_file,this_seq_file) == 0:
		print "####################################"
		print "ERROR: sequence "+seq_name+" not found in file "+seq_file
		print "Bye bye"
		quit()
	formatdb_nuc(this_seq_file)
	data=blastN_and_parseOutput(mono_file, this_seq_file)
	data_filtred=select_monomers_on_evalue(data, blast_maximal_evalue)
	data_filtred=select_monomers_on_score(data_filtred, blast_minimal_score)
	data_filtred=select_monomers_on_length(data_filtred, blast_minimal_length)
	data_no=get_nonoverlapping_monomers(data_filtred, blast_maximal_overlap)
	
	seq_monomers_data=insert_newhit_by_smithwaterman_in_empty_data(seq_monomers_data, ss, 1, data_seq[0]['len_seq'], seq_file, mono_file,waterman_minimal_score, waterman_minimal_length)







	
	
######################################

parser = argparse.ArgumentParser(description='Search for tandemly repeated monomers ...')
parser.add_argument('-s', dest="seq_file")
parser.add_argument('-m', dest="mono_file")
parser.add_argument('-o', dest="output_file")
parser.add_argument('-proc',  dest="processors",       type=int,         default=4)
parser.add_argument('-verb',  dest="verbose",                   default=10)
parser.add_argument('-bmo',  dest="blast_maximal_overlap",     type=int, default=5)
parser.add_argument('-bmd',  dest="max_dist_hits",             type=int, default=30)

parser.add_argument('-bml',  dest="blast_minimal_length",    type=int, default=100)
parser.add_argument('-bMe',  dest="blast_maximal_evalue",    type=float, default=10000.0)
parser.add_argument('-bms',  dest="blast_minimal_score",     type=int, default=0)


parser.add_argument('-wml',  dest="waterman_minimal_length", type=int, default=100)
parser.add_argument('-wms',  dest="waterman_minimal_score",  type=float, default=200)
parser.add_argument('-wgo',  dest="waterman_extension_gap_penality", type=float, default=10)
parser.add_argument('-wge',  dest="waterman_opening_gap_penality",   type=float, default=0.5)
args = parser.parse_args()


#############################
# creating blast databases ...

#############################

for record in SeqIO.parse(args.seq_file, "fasta"):
	sequence_analysis(record.id, args.seq_file, args.mono_file, args.blast_maximal_overlap, args.blast_minimal_length,args.blast_minimal_score, args.blast_maximal_evalue, args.waterman_minimal_score, args.waterman_minimal_length)

quit()

#############################
# geting sequences length ...
# reading blast output ...
print "reading data ..."
data = np.genfromtxt(args.blast_file+'4', delimiter=' ', usecols=[0,1,2,3,6,7,8,9,11,12,13,14],
      names=['seq','monomer','simil', 'length','mono_beg', 'mono_end', 'seq_beg', 'seq_end', 'score', 'len_mono', 'len_seq', 'strand']  , 
      dtype=None)
      #dtype=[('seq',object) ,('monomer',object), ('simil','f'), ('length','i') ,('mono_beg','i'), ('mono_end','i'), ('seq_beg','i'), ('seq_end','i'), ('score','f'), ('len_mono','i'), ('len_seq','i'), ('strand',object)])
#data=np.array(data, dtype=[('seq',object) ,('monomer',object), ('simil','f'), ('length','i') ,('mono_beg','i'), ('mono_end','i'), ('seq_beg','i'), ('seq_end','i'), ('score','f'), ('len_mono','i'), ('len_seq','i'), ('strand',object)])


workers = args.processors
work_queue = Queue()
done_queue = Queue()
processes = []

for ss in set(data['seq']):
	work_queue.put(ss)	

q = Process(target=writer_to_file, args=(done_queue,args.seq_file,args.output_file))
q.start()

for w in xrange(workers):
        p = Process(target=worker, args=(work_queue, done_queue, data, args.seq_file, args.mono_file, args.blast_maximal_overlap, args.blast_minimal_length, args.waterman_minimal_score, args.waterman_minimal_length))
        p.start()
        processes.append(p)
        work_queue.put('STOP')

for p in processes:
        p.join()

done_queue.put('STOP')
		
print "Bye bye !!!"

