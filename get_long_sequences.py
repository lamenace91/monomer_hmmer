#! /usr/bin/env python

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
parser.add_argument('-i', dest="seq_file")
parser.add_argument('-o', dest="output_file")
parser.add_argument('-t', dest="min_length",    type=int, default=50)
args = parser.parse_args()

if os.path.exists(args.output_file):
    os.remove(args.output_file)

seq_dict = SeqIO.index(args.seq_file, "fasta")

long_sequences = [] # Setup an empty list
handle=open(args.output_file, 'w')

for ss in seq_dict.keys():
	if len(seq_dict[ss]) >= args.min_length:
		long_sequences.append(seq_dict[ss])

print "input sequences: ", len(seq_dict)
print "output sequences: ", len(long_sequences)
handle=open(args.output_file, 'w')
SeqIO.write(long_sequences, handle, "fasta")
handle.close()
print "Bye bye !!!"
