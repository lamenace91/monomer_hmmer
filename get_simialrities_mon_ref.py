#! /usr/bin/env python

import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import my_function as mf

parser = argparse.ArgumentParser(description='Compute similarities between multiple sequences (monomers) and a reference sequence ...', epilog='')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-r', dest="ref_file", required=True)
parser.add_argument('-o', dest="out_file", required=True)
parser.add_argument('-a', dest="align_method", default="ggsearch")
parser.add_argument('-O', dest="options", default="")

args = parser.parse_args()

if not os.path.exists(args.seq_file):
	print("Error: sequence file ("+args.seq_file+") not found !!")
	quit()
if not os.path.exists(args.ref_file):
	print("Error: sequence file ("+args.ref_file+") not found !!")
	quit()
out_handle=open(args.out_file, "w")

handle=open(args.seq_file, "rU")
reference = SeqIO.read(open(args.ref_file), "fasta")
print ("seq1  ref id1 id2")
for record in SeqIO.parse(handle, "fasta") :
			id1, id2 = mf.get_similarity2(record.seq, reference.seq)
			print ("%s  %s %f %f" % (record.id, reference.id,  id1,id2))
out_handle.close()

