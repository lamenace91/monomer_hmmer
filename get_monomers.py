#! /usr/bin/python3.4


import os
import time
import argparse
import get_monomers_lib as gm



##################################################################################	
###################################################################################					
###################################################################################					
###################################################################################					
###################################################################################					
	
#~ def are_hits_overlapping(hit1, hit2):
	#~ if hit1['
###################################################################################					


parser = argparse.ArgumentParser(description='Search for tandem repeats ...', epilog='')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-r', dest="ref_file", required=True)
parser.add_argument('-o', dest="out_file", required=True)
parser.add_argument('-v',  dest="verbose",   type=int, default=0)

args = parser.parse_args()

if not os.path.exists(args.seq_file):
	print("Error: sequence file ("+args.seq_file+") not found !!")
	quit()
if not os.path.exists(args.ref_file):
	print("Error: sequence file ("+args.ref_file+") not found !!")
	quit()
	

#################################################################
if args.verbose > 0:
	print("Formating the hmmer db ...")
	t1=time.time()
monomer_hmm_file = gm.format_hmmdb(args.ref_file)
if args.verbose > 2:
	print("     hmm db in %s (%f s)" % (monomer_hmm_file, time.time()-t1))


#################################################################
if args.verbose > 0:
	print("Searching the sequences ...")
	t1=time.time()
hmmer_ouput_file = gm.search_hmm(args.seq_file, monomer_hmm_file)
if args.verbose > 2:
	print("     nhmmscan output in %s (%f s)" % (hmmer_ouput_file, time.time()-t1))

#################################################################
if args.verbose > 0:
	print("Parsing the hmmer output ...")
	t1=time.time()
hmm_hits = gm.parse_hmm_ouput(hmmer_ouput_file)
if args.verbose > 2:
	print("     number of hits: %d (%f s)" % (len(hmm_hits), time.time()-t1))
#~ for target_name in hmm_hits.keys():
	#~ alito=1
	#~ envto=1
	#~ for hit in hmm_hits[target_name]:
		#~ print("%s -> %d  %d   %d -- %d  %d -- %d (%s)" % (target_name, hit['alifrom']-alito-1,  hit['envfrom']-envto-1, hit['alifrom'], hit['alito'], hit['envfrom'], hit['envto'], hit["strand"]))
		#~ alito=hit['alito']
		#~ envto=hit['envto']
#print(hmm_hits)




