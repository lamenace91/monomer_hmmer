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

from Bio.SeqRecord import SeqRecord
import my_function as mf


##########################################
def writer_to_file(done_queue, seq_file, output_file, output_file_tab, output_file_rc, output_file_rc_seq, rc_index=1, verbose=0, naming=1):    
    print "Starting writing process ..."
    while 1:
        seq, m, rc = done_queue.get()
	if seq != 'STOP':
		handle=open(output_file_tab, 'a')
		print_analysis(m)
		if verbose > 5:
			print ("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
			print(m)
			print ("XXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
		print_ascii(seq,m)
		for ii in range(1,len(m)-1):
			handle.write(seq +" "+ m[ii][0] +" "+ str(m[ii][1]) +" "+ str(m[ii][2]) +" "+ str(m[ii][3]) +" "+   str(m[ii][4]) +" "+ str(m[ii][5]) +" "+ str(m[ii][6]) +" "+ str(m[ii][7]) +" "+   str(m[ii][8]) +" "+ str(m[ii][9])+" "+ str(m[ii][10]) + "\n")
			if naming == 2:
				mf.extract_sequence_to_newFile(seq,  seq_file, output_file, m[ii][5], m[ii][6], 'a', suffixe=m[ii][0])
			else:
				mf.extract_sequence_to_newFile(seq,  seq_file, output_file, m[ii][5], m[ii][6], 'a')	
		handle.close()
		
		# writing rolling_circle information
		if output_file_rc != "":
			handle=open(output_file_rc, 'a')
			handle.write(seq +" ")
			for cc in rc:
				handle.write( "," + str(int(cc)))
			handle.write("\n")
			handle.close()
			
			
			
		if output_file_rc_seq != "":		
			ss=mf.get_sequence(seq,seq_file).seq
			handle=open(output_file_rc_seq+str(rc_index), 'a')
			tmp=0
			ok=0
			for ii in range(0,len(ss)):
					if ok == 0 or rc[ii] > 0 and (rc[ii] == rc_index or (tmp < rc_index and rc[ii] >=rc_index) or (tmp > rc[ii] and rc[ii] >=rc_index)):
						if ok == 1:
							handle.write( "\n")
							print ""
						print ">"+seq+"_"+str(ii)
						handle.write( ">"+seq+"_"+str(ii)+"\n")
						ok=1
					handle.write((ss[ii]))
					print ss[ii],
					if rc[ii] > 0:
						tmp=rc[ii]
			handle.write("\n")
			print ""
			handle.close()
			
	else:
		break
    print "Ending writing process!!!"
    return True

######################################
def read_dict_of_sequences(seq_file, seq_format="fasta"):
	handle = open(seq_file, "rU")
	seq_dict = SeqIO.to_dict(SeqIO.parse(handle, seq_format))
	handle.close()
	return(seq_dict)


######################################
def read_rc_data(rc_file):
	rc={}
	handle = open(rc_file, "rU")
	for line in iter(handle):
		(seq,rc_line)=line.split()
		rc[seq]=map (int, rc_line.split(",")[1:])
	handle.close()
	return(rc)

	
######################################
def write_mono_from_rc_data(rc, seq_dict, output_file, rc_index=1):
	handle=open(output_file, 'a')

	for seq in rc.keys():
	#	print rc[seq]
		tmp=0
		ok=0
		for ii in range(0,len(seq_dict[seq])):
			if ok == 0 or rc[seq][ii] > 0 and (rc[seq][ii] == rc_index or (tmp < rc_index and rc[seq][ii] >=rc_index) or (tmp > rc[seq][ii] and rc[seq][ii] >=rc_index)):
				if ok == 1:
					handle.write( "\n")
					sys.stdout.write("\n")
				sys.stdout.write(">"+seq+"_"+str(ii)+"\n")
				handle.write( ">"+seq+"_"+str(ii)+"\n")
				ok=1
			handle.write((seq_dict[seq][ii]))
			sys.stdout.write(seq_dict[seq][ii])
			if rc[seq][ii] > 0:
				tmp=rc[seq][ii]
		handle.write("\n")
		print ""

	handle.close()

######################################
parser = argparse.ArgumentParser(description='Extract monomers by using rolling circle system ...', epilog='This program searchs for monomer limits into long sequences (-s option) by using the Smith and Waterman algorithm as implemented in _water_ (package EMBOSS, needed) and the monomers sequences (-m option)')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-r', dest="seq_rc_file", required=True)
parser.add_argument('-o', dest="output_file", required=True)
parser.add_argument('-i', dest="rc_index", type=int, default=1)
parser.add_argument('-verb',  dest="verbose",  type=int, default=10)
args = parser.parse_args()

if os.path.exists(args.output_file):
    os.remove(args.output_file)

#############################

## reading rolling circle data
rc=read_rc_data(args.seq_rc_file)

## reading sequences
seq_dict=read_dict_of_sequences(args.seq_file)

## write monomers
write_mono_from_rc_data(rc, seq_dict, args.output_file, args.rc_index)

print "Bye bye !!!"

