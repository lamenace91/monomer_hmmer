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
from Bio.Seq import Seq
import math
#from scipy import DataFrame
from multiprocessing import Lock, Process, Queue, current_process
from Bio.Align import AlignInfo


parser = argparse.ArgumentParser(description=' ...')
parser.add_argument('-i', dest="align_file")
parser.add_argument('-o', dest="output_file")
parser.add_argument('-n', dest="name", default='consensus')
parser.add_argument('-t', dest="threshold", type=float, default=0.7)
parser.add_argument('-c', dest="char",  default='N')

args = parser.parse_args()

if os.path.exists(args.output_file):
    os.remove(args.output_file)


f=open(args.align_file,"r")
alignment = AlignIO.read(f, "fasta")
f.close()
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus(args.threshold, ambiguous=args.char)
consensus.id='cons'
f=open(args.output_file,"w")
f.write(">"+args.name+"\n")
f.write(str(consensus)+'\n')
f.close()


