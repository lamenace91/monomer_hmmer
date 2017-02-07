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





def get_rolling_alignments(seq, mono, seq_file, mono_file, sw_ext, sw_open, verbose=0):
	#############################
	#### temporary_file
	align_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	mono_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	
	## getting sequence monomer
	mf.extract_sequence_to_newFile(seq.id,seq_file,seq_tmpfile, beg=mono[5], end=mono[6])
#	print(mono)
#	print(seq.id)
#	print(mono[5])
#	print(mono[6])

	# getting reference monomer	
	if mono[7] == 'plus':
		seqq=mf.get_sequence(mono[0],mono_file).seq
	else:
		seqq=mf.get_sequence(mono[0],mono_file).reverse_complement().seq	
	output_handle = open(mono_tmpfile, "w")
	SeqIO.write(SeqRecord(seqq,id="mono_tmp_name"), output_handle, "fasta")
	output_handle.close()
	
	# alignment
	os.system("needle  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 20 -gapextend 2 -aformat3 fasta 1>2 2> /dev/null ")

	# getting alignment
#	print(align_tmpfile)
#	print(mono_tmpfile)
#	print(seq_tmpfile)
	records=SeqIO.parse(open(align_tmpfile, "rU"), "fasta")
	
	mono_seq=records.next()
	seq_seq =records.next()
	#print seq.id
	#print mono[5]
	#print mono[6]
	#print mono_seq.seq
	#print seq_seq.seq
	rollingCircle=[]
	mm=0
	ss=0
	for ii in range(0, len(mono_seq.seq)):
		if mono_seq.seq[ii] != "-" and seq_seq.seq[ii] != "-" :
			mm=mm+1	
			ss=ss+1
#			print "yyy "+str(mm)+" "+seq_seq.seq[ii]+" "+str(ss)
			rollingCircle.append(mm)
		elif mono_seq.seq[ii] == "-" and seq_seq.seq[ii] != "-" :
			ss=ss+1
#			print "xxx "+str(0)+" "+seq_seq.seq[ii]+" "+str(ss)
			rollingCircle.append(0)
		elif mono_seq.seq[ii] != "-" and seq_seq.seq[ii] == "-" :
			mm=mm+1			
	#print align_tmpfile
	os.remove(align_tmpfile)
	os.remove(seq_tmpfile)
	os.remove(mono_tmpfile)
#	print rollingCircle
#	print mono[5]
	rollingCircle=np.array(rollingCircle)
	rollingCircle=rollingCircle
	
	
	return(rollingCircle)

################################
def re_estimation_monomer_limits(seq, mono1, mono2, seq_file, mono_file, sw_ext, sw_open):
	#############################
	#### temporary_file
	align_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	mono_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	if mono1[5] < mono2[6]:
		mf.extract_sequence_to_newFile(seq.id,seq_file,seq_tmpfile, beg=mono1[5], end=mono2[6])
	else:
		print "error: bla bla bla ... " + str(mono1[5]) + " " + str(mono2[6])
		print mono1
		print mono2
		quit()
		return(new_end, new_beg)
		
	seq=""
	if mono1[7] == 'plus':
		seq=mf.get_sequence(mono1[0],mono_file).seq
	else:
		seq=mf.get_sequence(mono1[0],mono_file).reverse_complement().seq
		#os.system("fastacmd -s "+mono1['monomer']+" -S 2 -d "+mono_file+">"+mono_tmpfile)
	if mono2[7] == 'plus':
		seq=seq+mf.get_sequence(mono2[0],mono_file).seq
		#os.system("fastacmd -s "+mono2['monomer']+" -d "+mono_file+" |grep -v \">\">>"+mono_tmpfile)
	else:
		zz=mf.get_sequence(mono2[0],mono_file)
		seq=seq+mf.get_sequence(mono2[0],mono_file).reverse_complement().seq
		#os.system("fastacmd -s "+mono2['monomer']+" -S 2 -d "+mono_file+" |grep -v \">\">>"+mono_tmpfile)
		
	output_handle = open(mono_tmpfile, "w")
	SeqIO.write(SeqRecord(seq,id="mono_tmp_name"), output_handle, "fasta")
	output_handle.close()
	
#	os.system("needle  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 2 -aformat3 fasta 2>1 > /dev/null")
#	os.system("needle  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 2 -aformat3 fasta &> /dev/null")
	os.system("needle  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 20 -gapextend 2 -aformat3 fasta 1>2 2> /dev/null ")
	#os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen "+str(sw_open)+ " -gapextend "+ str(sw_ext)+ " 1>2 2> /dev/null")

	f=open(align_tmpfile,"r")
	alignment = AlignIO.read(f, "fasta")
	f.close()
	os.remove(align_tmpfile)
	os.remove(seq_tmpfile)
	os.remove(mono_tmpfile)

	p1=0
	p2=0
	s1=0
	while p1 < mono1[9]:
		if alignment[0].seq[p2] != '-':
			p1+=1
		if alignment[1].seq[p2] != '-':
			s1+=1
		p2+=1
	new_end=mono1[5]+s1-1
	while p1 <= mono1[9]:
		if alignment[0].seq[p2] != '-':
			p1+=1
		if alignment[1].seq[p2] != '-':
			s1+=1
		p2+=1
	if alignment[1].seq[p2-1] == '-':
		s1+=1
	new_beg=mono1[5]+s1-1
	
	
	return(new_end, new_beg)


######################################
## get new hit from smith waterman alignement
######################################
def newhit_by_smithwaterman(seq, beg, end, mono_file, sw_ext, sw_open,verbose=0):
#	print "         searching new hit by SW ..."
	new_mono_monomer=''
	score_max=0.0
	new_mono_simil=-1
	new_mono_mono_beg=-1
	new_mono_mono_end=-1
	new_mono_seq_beg=-1
	new_mono_seq_end=-1
	new_mono_length=-1
	new_mono_score=-1
	new_mono_monomer=''
	new_mono_strand=''
	score_max=-1
	new_mono_seq=''
	new_mono_mono=''
	new_mono_length=0	
	align_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	mono_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
#	print "tmp files are:"
#	print align_tmpfile
#	print seq_tmpfile
#	print mono_tmpfile
	if beg < end:
#		print "          writing sequence ..."
		output_handle = open(seq_tmpfile, "w")
		SeqIO.write(SeqRecord(seq.seq[beg:end],id=seq.id), output_handle, "fasta")
		output_handle.close()
#		print "                            done"
	else:
		print "begin has to be lower than end !!"
		return [new_mono_monomer,new_mono_simil,new_mono_length,new_mono_mono_beg,new_mono_mono_end,beg+new_mono_seq_beg-1,beg+new_mono_seq_end-1,new_mono_strand,new_mono_score, new_mono_length]
#	print "          reading monomers ..."
	new_found=0
	handle = open(mono_file, "rU")
	record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
#	print "                            done"
		
#	print "          looping on monomers ..."
	for mono_record in record_dict:
#		print "              writing monomer in tmp file: "+mono_tmpfile
		output_handle = open(mono_tmpfile, "w")
		SeqIO.write(SeqRecord(record_dict[mono_record].seq,id=mono_record), output_handle, "fasta")
		output_handle.close()
#		print "                            done"
#		print "              aligning... "
#		os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen "+str(sw_open)+ " -gapextend "+ str(sw_ext)+ " 2>1 > /dev/null")
#		os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen "+str(sw_open)+ " -gapextend "+ str(sw_ext)+ " &> /dev/null")
		os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen "+str(sw_open)+ " -gapextend "+ str(sw_ext)+ " 1>2 2> /dev/null")
#		print "                        done"
#		print "              parsing... "
		(mono_seq, seq_seq, mono_beg,mono_end,seq_beg,seq_end, align_length,score)=mf.parse_water_alignment(align_tmpfile,verbose)		
#		print (mono_seq, seq_seq, mono_beg,mono_end,seq_beg,seq_end, align_length,score)
#		print "                        done"
		if score >= score_max:
#			print "                NEW MONO !!!"
			new_found=1
			new_mono_simil=round(mf.get_similarity(mono_seq,seq_seq)*100,3)
			new_mono_mono_beg=mono_beg
			new_mono_mono_end=mono_end
			new_mono_seq_beg=seq_beg
			new_mono_seq_end=seq_end
			new_mono_length=align_length
			new_mono_score=score
			new_mono_monomer=record_dict[mono_record].id
			new_mono_strand='plus'
			new_mono_length=align_length
			new_mono_seq=seq_seq
			new_mono_mono=mono_seq
			score_max=score
			new_mono_length=len(record_dict[mono_record].seq)
#		print "              writing REVERSE monomer in tmp file: "+mono_tmpfile
		output_handle = open(mono_tmpfile, "w")
		SeqIO.write(SeqRecord(record_dict[mono_record].seq.reverse_complement(),id=mono_record), output_handle, "fasta")
		output_handle.close()
#		print "                            done"
#		print "              aligning... "		
#		os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen "+str(sw_open)+ " -gapextend "+ str(sw_ext)+ " 2>1 > /dev/null")
#		os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen "+str(sw_open)+ " -gapextend "+ str(sw_ext)+ " &> /dev/null")
		os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen "+str(sw_open)+ " -gapextend "+ str(sw_ext)+ " 1>2 2> /dev/null")
#		print "                        done"
#		print "              parsing... "
		(mono_seq, seq_seq, mono_beg,mono_end,seq_beg,seq_end, align_length,score)=mf.parse_water_alignment(align_tmpfile,verbose)
#		print (mono_seq, seq_seq, mono_beg,mono_end,seq_beg,seq_end, align_length,score)
#		print "                        done"

#		simil=mf.get_similarity(mono_seq,seq_seq)
		if score >= score_max:
			#print "                NEW MONO !!!"
			new_found=1
			new_mono_simil=round(mf.get_similarity(mono_seq,seq_seq)*100,3)
			new_mono_mono_beg=mono_beg
			new_mono_mono_end=mono_end
			new_mono_seq_beg=seq_beg
			new_mono_seq_end=seq_end
			new_mono_length=align_length
			new_mono_score=score
			new_mono_monomer=record_dict[mono_record].id
			new_mono_strand='minus'
			score_max=score
			new_mono_seq=seq_seq
			new_mono_mono=mono_seq
			new_mono_length=len(record_dict[mono_record].seq)
	os.remove(align_tmpfile)
	os.remove(seq_tmpfile)
	os.remove(mono_tmpfile)

	return [new_mono_monomer,new_mono_simil,new_mono_length,new_mono_mono_beg,new_mono_mono_end,beg+new_mono_seq_beg-1,beg+new_mono_seq_end-1,new_mono_strand,new_mono_score,new_mono_length]
	


######################################
## Run analysis for one sequence
######################################
def sequence_analysis(work_queue, done_queue, seq_file, mono_file, output_file_rc,
							waterman_minimal_score, waterman_minimal_length, waterman_minimal_simil, 
							sw_ext, sw_open,verbose=0):
    print "Starting analysis process "+current_process().name+" ..."
    for seq_name in iter(work_queue.get, 'STOP'):
	print "#########################################"
	print current_process().name+" : "+seq_name
        sys.stdout.flush()
	seq=mf.get_sequence(seq_name,seq_file)
	seqlen=len(seq.seq)
	seqid=seq.id
	seqid.split("_")
#	print "Inserting fake hits ... (0,end)"
	hits=[]
	mono_beg=['monomer', 0.0, 0, 0, 0, 0, 0,'na', 0.0, 0, 0]
	mono_end=['monomer', 0.0, 0, 0, 0, seqlen, seqlen,'na', 0.0, 0, 1]
	hits.append(mono_beg)
	hits.append(mono_end)
	ii=1
#	print "Looping between hits (n="+str(len(hits))+")"
	
	while ii < len(hits):
		#print "	   looking at hits: "+str(ii-1)+" "+str(ii)
		#print "	        located at: "+str(hits[ii-1][6])+" "+str(hits[ii][5])
		
		#print(hits)	
		if (hits[ii][5] - hits[ii-1][6]) > 10 and hits[ii][10] == 1:
	#		print "	        aligning..."
			if (hits[ii][5] - hits[ii-1][6]) < 1000:
				new_hit=newhit_by_smithwaterman(seq, hits[ii-1][6], hits[ii][5], mono_file, sw_ext, sw_open,verbose)
			else:
				seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
				mf.extract_sequence_to_newFile(seqid,seq_file, seq_tmpfile, hits[ii-1][6], hits[ii][5])
				mf.formatdb_nuc(seq_tmpfile)
				blast_outfile=mf.blast_m8(mono_file,seq_tmpfile)
				blast_data=np.atleast_1d(mf.read_blast(blast_outfile))
				os.remove(seq_tmpfile)
				os.remove(blast_outfile)
			#	print(blast_outfile)
				#print(blast_data)
				#exit()
				if len(blast_data) < 1:
					new_hit=newhit_by_smithwaterman(seq, 0, -1, mono_file, sw_ext, sw_open,verbose)
				else:
					if blast_data[0][8] <  blast_data[0][9]:
						beg=max(hits[ii-1][6]+blast_data[0][8]-200, hits[ii-1][6])
						end=min(hits[ii-1][6]+blast_data[0][9]+200, hits[ii][5])
					else:
						beg=max(hits[ii-1][6]+blast_data[0][9]-200, hits[ii-1][6])
						end=min(hits[ii-1][6]+blast_data[0][8]+200, hits[ii][5])
					#print ("toto")
					#print(beg)
					#print(end)
				
					new_hit=newhit_by_smithwaterman(seq, beg, end, mono_file, sw_ext, sw_open,verbose)
				
		#	print "	                 done!!!!!!"
			if verbose > 7:
				print ("------------------------------------------")
				print (new_hit)
				print waterman_minimal_length
				print waterman_minimal_score
				print waterman_minimal_simil
				print ("------------------------------------------")
				
			if new_hit[8] > waterman_minimal_score and new_hit[2] > waterman_minimal_length and new_hit[1] > waterman_minimal_simil:
				if verbose > 7:
					print("added")
				new_hit.append(1)
				hits.append(new_hit)
				hits.sort(key=lambda x: x[5])
				ii=ii-1
			else:
				if verbose > 7:
					print("not added")
				
		ii=ii+1		
#	print (hits)	
		
#	print "//////////////////////////////"
#	print "//////////////////////////////"
#	print "//////////////////////////////"
#	print "//////////////////////////////"
#	print "Looping on all previous hits"
#	print "//////////////////////////////"
#	print "//////////////////////////////"
#	print "//////////////////////////////"
#	print "//////////////////////////////"
		
	for ii in range(1,len(hits)-2):
		
		dist_seq = hits[ii+1][5]-hits[ii][6]
		dist_mono=(hits[ii][9] -hits[ii][4])+hits[ii+1][3]
#		print dist_seq
#		print dist_mono
       		sys.stdout.flush()
		
		
		if dist_seq != 1 and (math.fabs(dist_mono-dist_seq) < 30 or dist_seq < 30):
			#print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
			#print(hits[ii]) 
			#print(hits[ii+1])
			#print("~~~~")
			new_end, new_beg=re_estimation_monomer_limits(seq, hits[ii], hits[ii+1], seq_file, mono_file, sw_ext, sw_open)
			hits[ii][6]=new_end
			hits[ii+1][5]=new_beg
	
			#print hits[ii]
			#print hits[ii+1]
			#print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
##
	#print "|||||||||||||||||||||||||||||||||||||||"
	#print (hits)	
	ii=1
	L=len(hits)
	while ii < L:
		if hits[ii][5] >= hits[ii][6]:
			#print("del")
			#print(hits[ii])
			del hits[ii]
			ii=ii-1
		L=len(hits)
		ii=ii+1


	#print "|||||||||||||||||||||||||||||||||||||||"
	print (hits)	
	rc=np.array([])

	if output_file_rc != "":
		rc=np.array([])
		tmp=0
		for ii in range(1,len(hits)-1):
		#	print tmp, hits[ii][5], hits[ii][6]
			for jj in range(tmp, hits[ii][5]):
				rc=np.concatenate((rc, [-1]))
			gr= get_rolling_alignments(seq, hits[ii], seq_file, mono_file, sw_ext, sw_open, verbose)
			#print(hits[ii][6]- hits[ii][5]+1)
			#print(len(gr))
			#print(gr)
			rc=np.concatenate((rc, gr))
			tmp=hits[ii][6]+1
		for jj in range(tmp, seqlen):
			rc=np.concatenate((rc, [-2]))
	done_queue.put((seq_name, hits, rc))
    print "Ending analysis process "+current_process().name+" !!!"
    quit()
    return True























##########################################
def writer_to_file(done_queue, seq_file, output_file, output_file_tab, output_file_rc, output_file_rc_seq, rc_index=1, verbose=0, naming=1):    
    print "Starting writing process ..."
    while 1:
        seq, m, rc = done_queue.get()
	if seq != 'STOP':
		handle=open(output_file_tab, 'a')
		print(m)
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
		
		# wrinting rollong_circle information
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
	
##########################################
def print_ascii(seq, m): 
	print(m)
	m=np.atleast_2d(m)
	init_gap= m[1][5]
	print(m)

	output=""
	if init_gap > 0:
		output=output+ ".."+str(init_gap)+".."
	for ii in range(1,len(m)-1):
		leng=m[ii][6]-m[ii][5]
		if m[ii][7] == "plus":
			output=output+ "--"+str( leng)+"-->" 
		else:
			output=output+ "<-"+str( leng)+"--" 
		dist=m[ii+1][5]-m[ii][6]	
		if dist != 1 :
			output=output+".."+str(dist)+".."
	print output
			
	for ii in range(1,len(m)-1):
			print ' ',	
######################################
def print_analysis(data):
	print(data)
	print("######################################")
	total_gap=0
	#dt = np.dtype([np.str_, np.float, np.int, np.int, np.int, np.int, np.int, np.str_, np.float, np.int, np.int])
	
	dt= np.dtype([('monomer', 'a100'), ('score', 'f8'), ('l', 'i4'), ('b', 'i4'), ('c', 'i4'), ('d', 'i4'), ('e', 'i4'),('strand', 'a100'),('h', 'f8')  ,('f', 'i4'), ('g', 'i4') ])
	xx=np.asarray(data)
	print("######################################")
	print("######################################")
	print("######################################")
	print("######################################")
	print(xx)
	print("######################################")
	data=np.atleast_2d(xx)
	data.dtype=dt
	nb_lines=len(data)-1
	seq_len=int(data[nb_lines][5])
	print("--------------------------------------")
	print(data)
	print("######################################")
	print (seq_len)
	print(nb_lines)
	print("------------")
	print(seq_len-1)
	print(nb_lines-1)
	print("------------")
	init_gap= data[1][5]
	term_gap= seq_len-int(data[nb_lines-1][6])
	
	for ii in range(1,nb_lines):
		total_gap+=(data[ii][5]-data[ii-1][6]-1)
	total_gap+=term_gap
	print("seq len: "+str(seq_len)+" initial_gap: "+str(init_gap)+" term gap: "+str(term_gap)+" total gap: "+str(total_gap))
	return(data)


######################################
parser = argparse.ArgumentParser(description='Search for tandemly repeated monomers ...', epilog='This program searchs for monomer limits into long sequences (-s option) by using the Smith and Waterman algorithm as implemented in _water_ (package EMBOSS, needed) and the monomers sequences (-m option)')
parser.add_argument('-s', dest="seq_file", required=True)
parser.add_argument('-m', dest="mono_file", required=True)
parser.add_argument('-o', dest="output_file", required=True)
parser.add_argument('-O', dest="output_file_tab", required=True)
parser.add_argument('-r', dest="output_file_rc", default="")
parser.add_argument('-R', dest="output_file_rc_seq", default="")
parser.add_argument('-i', dest="rc_index", type=int, default=1)
parser.add_argument('-n', dest="naming",     type=int,         default=2)
parser.add_argument('-proc',  dest="processors",       type=int,         default=4)
parser.add_argument('-verb',  dest="verbose",          type=int,              default=10)
parser.add_argument('-wml',  dest="waterman_minimal_length", type=int, default=0)
parser.add_argument('-wms',  dest="waterman_minimal_score",  type=float, default=100)
parser.add_argument('-wmS',  dest="waterman_minimal_simil",  type=float, default=0.0)
parser.add_argument('-wge',  dest="waterman_extension_gap_penality", type=float, default=1)
parser.add_argument('-wgo',  dest="waterman_opening_gap_penality",   type=float, default=10)




args = parser.parse_args()

if os.path.exists(args.output_file):
    os.remove(args.output_file)
if os.path.exists(args.output_file_tab):
    os.remove(args.output_file_tab)
if args.output_file_rc != "" and os.path.exists(args.output_file_rc):
    os.remove(args.output_file_rc)
if args.output_file_rc_seq != "" and os.path.exists(args.output_file_rc_seq+str(args.rc_index)):
    os.remove(args.output_file_rc_seq+str(args.rc_index))


#############################
mf.formatdb_nuc( args.mono_file)

workers = args.processors
work_queue = Queue()
done_queue = Queue()
processes = []

N=0
for record in SeqIO.parse(args.seq_file, "fasta"):
	N=N+1
	work_queue.put(record.id)	
if args.verbose > 2:
	print "number of sequences: "+str(N)
q = Process(target=writer_to_file, args=(done_queue,args.seq_file,args.output_file,args.output_file_tab, args.output_file_rc, args.output_file_rc_seq,  args.rc_index, args.verbose, args.naming))
q.start()

for w in xrange(workers):
        p = Process(target=sequence_analysis, args=(work_queue, done_queue, args.seq_file, args.mono_file, args.output_file_rc, args.waterman_minimal_score, args.waterman_minimal_length, args.waterman_minimal_simil,args.waterman_extension_gap_penality, args.waterman_opening_gap_penality,args.verbose))						
        p.start()
        processes.append(p)
        work_queue.put('STOP')

for p in processes:
        p.join()


done_queue.put(('STOP','STOP','STOP'))

print "Waiting for the writing  process ..."
q.join()
	
print "Bye bye !!!"

