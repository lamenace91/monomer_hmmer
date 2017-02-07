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

#import threading
#############################################################



def insert_newhit_by_smithwaterman_in_empty_data(data, seq, beg, end, seq_file, mono_file, waterman_minimal_score, waterman_minimal_length):
	(monomer,simil,length,mono_beg,mono_end,seq_beg,seq_end,strand,score)=newhit_by_smithwaterman(seq, beg, end, args.seq_file, args.mono_file)
	if monomer == '':
		return data
	if monomer != "" and score > waterman_minimal_score and length > waterman_minimal_length:
		data=np.append(data[0],data[0])
		data=np.delete(data,1,axis=0)
		last=len(data)-1
		data[last]['monomer']=monomer
		data[last]['simil']=simil
		data[last]['length']=length
		data[last]['seq_beg']=seq_beg
		data[last]['seq_end']=seq_end
		data[last]['mono_beg']=mono_beg
		data[last]['mono_end']=mono_end
		data[last]['score']=score
		data.sort(order='seq_beg')
	return data










#############################################################
def read_array(filename, dtype, separator=',', col=None):
       """ Read a file with an arbitrary number of columns.
           The type of data in each column is arbitrary
           It will be cast to the given dtype at runtime
       """
       cast = np.cast
       data = [[] for dummy in xrange(len(dtype))]
       if col == None:
       		col=xrange(len(dtype))
       print col
       for line in open(filename, 'r'):
           
	  fields = line.strip().split(separator)
	  nb=0
          for i, number in enumerate(fields):
	      if i in col:
              	 data[nb].append(number)
		 nb=nb+1
       print "####################"
       for i in xrange(len(dtype)):
           data[i] = cast[dtype[i]](data[i])
       print data
       return np.rec.array(data, dtype=dtype)
       
#############################################################
def num(seq, beg, end, seq_file, mono_file, waterman_minimal_score, waterman_minimal_length):
	print "okkkk"
	data=[]
	(monomer,simil,length,mono_beg,mono_end,seq_beg,seq_end,strand,score)=newhit_by_smithwaterman(seq, beg, end, args.seq_file, args.mono_file)
	if monomer == '':
		return data
	print monomer
	print score
	print length
	print waterman_minimal_score
	print waterman_minimal_length
	if monomer != "" and score > waterman_minimal_score and length > waterman_minimal_length:
#		data=np.append(data[0],data[0])
#		data=np.delete(data,1,axis=0)
#		last=len(data)-1
		data.append({})
		last=0
		data[last]['monomer']=monomer
		data[last]['simil']=simil
		data[last]['length']=length
		data[last]['seq_beg']=seq_beg
		data[last]['seq_end']=seq_end
		data[last]['mono_beg']=mono_beg
		data[last]['mono_end']=mono_end
		data[last]['score']=score
		print data
		data.sort(order='seq_beg')
	print data
	return data


#############################################################
def insert_newhit_by_smithwaterman(data, ii, jj, seq_file, mono_file, waterman_minimal_score, waterman_minimal_length):
	if ii >=0 and jj < 0:
	
		print "terminal ...."
		seq=data[ii]['seq']
		beg=data[ii]['seq_end']
		end=data[ii]['len_seq']
	elif ii < 0 and jj >= 0:
		print "initial ...."
		seq=data[ii]['seq']
		beg=1
		end=data[jj]['seq_beg']
	elif ii >= 0 and jj >= 0:
		print "internal ...."
		seq=data[ii]['seq']
		beg=data[ii]['seq_end']
		end=data[jj]['seq_beg']

	(monomer,simil,length,mono_beg,mono_end,seq_beg,seq_end,strand,score)=newhit_by_smithwaterman(seq, beg, end, args.seq_file, args.mono_file)
	if monomer == '':
		return data
	if monomer != "" and score > waterman_minimal_score and score > waterman_minimal_length:
		data=np.append(data,data[0])
		last=len(data)-1
		data[last]['seq']=seq
		data[last]['monomer']=monomer
		data[last]['simil']=simil
		data[last]['length']=length
		data[last]['seq_beg']=seq_beg
		data[last]['seq_end']=seq_end
		data[last]['mono_beg']=mono_beg
		data[last]['mono_end']=mono_end
		data[last]['score']=score
		data[last]['strand']=strand
		data.sort(order='seq_beg')
	return data
#############################################################
def parse_water_alignment(alignfile):
	score=0.0
	f=open(alignfile,'r')
	lines=f.xreadlines()
#	lines = [line.strip() for line in open(alignfile)]
#	lines = [l for l in f.readlines() if l.strip()]

	lines = filter(lambda x: not re.match(r'^\s*$', x),lines )
	for ll in lines:
		ll=ll.strip()
		tmp=ll.split()
		if len(tmp)==3 and tmp[1]=="Score:":
			score=float(tmp[2])
	lines = filter(lambda x: not re.match(r'^#.*$', x),lines )
	ii=0
	seq1=''
	beg1=0
	end1=0
	seq2=''
	beg2=0
	end2=0
	for ll in lines:
		ll=ll.strip()
		tmp=ll.split()
		if (ii%3)==0:
			seq1=seq1+tmp[2]
			end1=tmp[3]
			if ii == 0:
				beg1=tmp[1]
		elif (ii%3)==2:
			seq2=seq2+tmp[2]
			end2=tmp[3]
			if ii == 2:
				beg2=tmp[1]
		ii+=1
	return(seq1, seq2, int(beg1), int(end1), int(beg2), int(end2), len(seq1),score)
#############################################################
def newhit_by_smithwaterman(seq, beg, end, seq_file, mono_file):
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
	align_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	mono_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	if beg < end:
		os.system("fastacmd -s "+seq+" -L"+str(beg)+","+str(end)+" -d "+seq_file+">"+seq_tmpfile)
	else:
		return (new_mono_monomer,new_mono_simil,new_mono_length,new_mono_mono_beg,new_mono_mono_end,beg+new_mono_seq_beg-1,beg+new_mono_seq_end-1,new_mono_strand,new_mono_score)
	new_found=0
	for mono_record in SeqIO.parse(mono_file, "fasta"):
		os.system("fastacmd -s "+mono_record.id+" -d "+mono_file+">"+mono_tmpfile)
		os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 10 2>1 > /dev/null")
		(mono_seq, seq_seq, mono_beg,mono_end,seq_beg,seq_end, align_length,score)=parse_water_alignment(align_tmpfile)
		if score >= score_max:
			new_found=1
			new_mono_simil=round(get_similarity(mono_seq,seq_seq)*100,3)
			new_mono_mono_beg=mono_beg
			new_mono_mono_end=mono_end
			new_mono_seq_beg=seq_beg
			new_mono_seq_end=seq_end
			new_mono_length=align_length
			new_mono_score=score
			new_mono_monomer=mono_record.id
			new_mono_strand='plus'
			new_mono_length=align_length
			new_mono_seq=seq_seq
			new_mono_mono=mono_seq
			score_max=score

		os.system("fastacmd -s "+mono_record.id+" -S 2 -d "+mono_file+">"+mono_tmpfile)
		os.system("water  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 10 2>1 > /dev/null")
		(mono_seq, seq_seq, mono_beg,mono_end,seq_beg,seq_end, align_length,score)=parse_water_alignment(align_tmpfile)
		simil=get_similarity(mono_seq,seq_seq)
		if score >= score_max:
			new_found=1
			new_mono_simil=round(get_similarity(mono_seq,seq_seq)*100,3)
			new_mono_mono_beg=mono_beg
			new_mono_mono_end=mono_end
			new_mono_seq_beg=seq_beg
			new_mono_seq_end=seq_end
			new_mono_length=align_length
			new_mono_score=score
			new_mono_monomer=mono_record.id
			new_mono_strand='minus'
			score_max=score
			new_mono_seq=seq_seq
			new_mono_mono=mono_seq

	os.remove(align_tmpfile)
	os.remove(seq_tmpfile)
	os.remove(mono_tmpfile)
	return (new_mono_monomer,new_mono_simil,new_mono_length,new_mono_mono_beg,new_mono_mono_end,beg+new_mono_seq_beg-1,beg+new_mono_seq_end-1,new_mono_strand,new_mono_score)

#############################################################
def get_limits_on_mono(seq1,seq2):
	ii=0
	beg=0
	while ii < len(seq2) and seq2[ii] == '-':
		if seq1[beg] != '-':
			beg+=1
		ii+=1
	pp=beg
	while ii < len(seq2):
		if seq2[ii] != '-':
			pp=ii
		ii+=1
	ii=0
	end=0
	while ii <= pp:
		if seq1[ii] != '-':
			end+=1
		ii+=1
	return (beg, end)
#############################################################
def get_similarity(seq1, seq2):
	simil=0.0
	total=0.0
	ii=0
	while ii < len(seq1):
		if seq1[ii] != '-' and seq2[ii] != '-':
			if seq1[ii] == seq2[ii]:
				simil+=1.0
			total+=1.0
		ii+=1
	if total > 0:
		return simil/total
	else:
		return -1.0

#############################################################
def align_mono(mono, seq_file, mono_file):
	new_mono=mono
	simil_max=0
	align_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	mono_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	if mono['seq_beg'] < mono['seq_end']:
		os.system("fastacmd -s "+mono['seq']+" -L"+str(mono['seq_beg'])+","+str(mono['seq_end'])+" -d "+seq_file+">"+seq_tmpfile)
	else:
		print "error: blu blu blu ..."
		quit()
	for mono_record in SeqIO.parse(mono_file, "fasta"):
   	 	print mono_record.id
		if mono['strand'] == 'plus':
			os.system("fastacmd -s "+mono_record.id+" -d "+mono_file+">"+mono_tmpfile)
		else:
			os.system("fastacmd -s "+mono_record.id+" -S 2 -d "+mono_file+">"+mono_tmpfile)
		os.system("needle  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 0.5 -aformat3 fasta 2>1 > /dev/null")
		f=open(align_tmpfile,"r")
          	try:
			alignment = AlignIO.read(align_tmpfile, "fasta")
  		except:
			print "#################"
			print "Error: .........."
			print mono
			print "file: "+align_tmpfile
			print "#################"
			return new_mono
		f.close()
		simil=get_similarity(alignment[0].seq,alignment[1].seq)
		if simil > simil_max:
			beg,end=get_limits_on_mono(alignment[0].seq,alignment[1].seq)
			new_mono['simil']=round(simil*100,3)
			new_mono['mono_beg']=beg
			new_mono['mono_end']=end
			new_mono['length']=len(alignment[0].seq)
			new_mono['monomer']=mono_record.id
	new_mono['score']=0.0
	os.remove(align_tmpfile)
	os.remove(seq_tmpfile)
	os.remove(mono_tmpfile)
	return new_mono

#############################################################
def 	get_monomers_from_small_hits(dat, seq, beg, end):
		new_mono=[]
		sub_data_tmp1=dat[dat[:]['seq_beg'] < data[:]['seq_end']]
		sub_data_tmp2=sub_data_tmp1[sub_data1[:]['seq']==seq]
		sub_data_tmp3=sub_data_tmp2[sub_data_tmp2[:]['seq_beg'] >=beg]
		sub_data_tmp4=sub_data_tmp3[sub_data_tmp3[:]['seq_end'] <=end]
		if len(sub_data_tmp4) == 0:
			return new_mono
		ind = np.argmax(sub_data_tmp4[:]['score'])
		new_mono2=sub_data_tmp4[ind]
		return new_mono2


#############################################################
def 	build_monomers_from_small_hits2(dat, seq, beg, end, overlap):
		new_mono=[]
		new_mono2=new_mono

		sub_data_tmp1=dat[dat[:]['seq_beg'] < data[:]['seq_end']]
		sub_data_tmp2=sub_data_tmp1[sub_data1[:]['seq']==seq]
		sub_data_tmp3=sub_data_tmp2[sub_data_tmp2[:]['seq_beg'] >=beg]
		sub_data_tmp4=sub_data_tmp3[sub_data_tmp3[:]['seq_end'] <=end]
		new_len=0
		if len(sub_data_tmp4) == 0:
			return new_mono
		for mm in set(sub_data_tmp4['monomer']):
			sub_data_tmp5=sub_data_tmp4[sub_data_tmp4[:]['monomer']==mm]
			for ss in set(sub_data_tmp5['strand']):
				sub_data_tmp=sub_data_tmp5[sub_data_tmp5[:]['strand']==ss]
				ind = np.argmax(sub_data_tmp[:]['score'])
				new_mono=sub_data_tmp[ind]
				sub_data_tmp[ind]['score']=-1
				selected_monomers=1
				new_mono['score']=1

				while selected_monomers < len(sub_data_tmp) and sub_data_tmp[ind]['score'] > 0:
					ind = np.argmax(sub_data_tmp[:]['score'])
					if new_mono['seq_end'] < (sub_data_tmp[ind]['seq_beg']+overlap) or new_mono['seq_beg'] > (sub_data_tmp[ind]['seq_end']-overlap):
						if new_mono['mono_end'] < (sub_data_tmp[ind]['mono_beg']+overlap) or new_mono['mono_beg'] > (sub_data_tmp[ind]['mono_end']-overlap):
							new_mono['mono_end']=max(new_mono['mono_end'], sub_data_tmp[ind]['mono_end'])
							new_mono['mono_beg']=min(new_mono['mono_beg'], sub_data_tmp[ind]['mono_beg'])
							new_mono['seq_end']=max(new_mono['seq_end'], sub_data_tmp[ind]['seq_end'])
							new_mono['seq_beg']=min(new_mono['seq_beg'], sub_data_tmp[ind]['seq_beg'])
							new_mono['score']+=1
							new_mono['simil']=(new_mono['simil']*new_mono['length']+sub_data_tmp[ind]['simil']*sub_data_tmp[ind]['length'])/(new_mono['length']+sub_data_tmp[ind]['length'])
							new_mono['length']=new_mono['length']+sub_data_tmp[ind]['length']
							sub_data_tmp[ind]['score']=-1
							selected_monomers=selected_monomers+1
#							print sub_data_tmp[ind]
						else:
							sub_data_tmp[ind]['score']=-1
					else:
						sub_data_tmp[ind]['score']=-1

				if new_mono['seq_end']-new_mono['seq_beg'] > new_len:
					new_mono2=new_mono
					new_len=new_mono['seq_end']-new_mono['seq_beg']
		return new_mono2



#############################################################
def 	build_monomers_from_small_hits(dat, seq, beg, end, overlap):
		new_mono=[]
		sub_data_tmp1=dat[dat[:]['seq_beg'] < data[:]['seq_end']]
		sub_data_tmp2=sub_data_tmp1[sub_data1[:]['seq']==seq]
		sub_data_tmp3=sub_data_tmp2[sub_data_tmp2[:]['seq_beg'] >=region_beg]
		sub_data_tmp=sub_data_tmp3[sub_data_tmp3[:]['seq_end'] <=region_end]
		if len(sub_data_tmp) == 0:
			return new_mono
		ind = np.argmax(sub_data_tmp[:]['score'])
		new_mono=sub_data_tmp[ind]
		sub_data_tmp[ind]['score']=-1
		sub_data_tmp=sub_data_tmp[sub_data_tmp[:]['monomer']==sub_data_tmp[ind]['monomer'] ]
		ind = np.argmax(sub_data_tmp[:]['score'])
		selected_monomers=1
		new_mono['score']=1
		while selected_monomers < len(sub_data_tmp) and sub_data_tmp[ind]['score'] > 0:
			ind = np.argmax(sub_data_tmp[:]['score'])
			if new_mono['seq_end'] < (sub_data_tmp[ind]['seq_beg']+overlap) or new_mono['seq_beg'] > (sub_data_tmp[ind]['seq_end']-overlap):
				if new_mono['mono_end'] < (sub_data_tmp[ind]['mono_beg']+overlap) or new_mono['mono_beg'] > (sub_data_tmp[ind]['mono_end']-overlap):
					new_mono['mono_end']=max(new_mono['mono_end'], sub_data_tmp[ind]['mono_end'])
					new_mono['mono_beg']=min(new_mono['mono_beg'], sub_data_tmp[ind]['mono_beg'])
					new_mono['seq_end']=max(new_mono['seq_end'], sub_data_tmp[ind]['seq_end'])
					new_mono['seq_beg']=min(new_mono['seq_beg'], sub_data_tmp[ind]['seq_beg'])
					new_mono['score']+=1
					new_mono['simil']=(new_mono['simil']*new_mono['length']+sub_data_tmp[ind]['simil']*sub_data_tmp[ind]['length'])/(new_mono['length']+sub_data_tmp[ind]['length'])
					new_mono['length']=new_mono['length']+sub_data_tmp[ind]['length']
					sub_data_tmp[ind]['score']=-1
					print sub_data_tmp[ind]
				else:
					sub_data_tmp[ind]['score']=-1
			else:
				sub_data_tmp[ind]['score']=-1
		return new_mono


#############################################################

def get_sequence(data, seq_file, output_file):
	if data['seq_beg'] < data['seq_end']:
		os.system("fastacmd -s "+data['seq']+" -L "+str(data['seq_beg'])+","+str(data['seq_end'])+" -S 1 -d "+seq_file+" >> "+ output_file)
	else:
		print "WARNING: "+ data['seq']+ " : " + str(data['seq_beg'])+"-"+str(data['seq_end'])
		return -1
	return 0

#############################################################
def get_nonoverlapping_monomers(data, blast_maximal_overlap, blast_minimal_length):
	nb_selected_monomers = 0
	selected_monomers=[]
	print "--------------------------------"
	print data
	data=data[data['length'] >= blast_minimal_length]
	print data
	print "okay        "
	if len(data) == 1:
		return data
	print "--------------------------------"	
	print data
	ind = np.argmax(data[0]['score'])
	print ind
	print "--------------------------------"	
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

	print selected_monomers_data
	
	if len(selected_monomers_data) > 1:
		selected_monomers_data.sort(order='seq_beg')
		
	return(selected_monomers_data)


################################
def check_monomers(data):
	ii=0
	while ii < len(data):
		if data[ii]['seq_beg'] > data[ii]['seq_end'] and data[ii]['seq'] != 'seq': 
			print data
			data=np.delete(data,ii,0)
			print "###################################"
			print "PROBLEM: monomer inverted"
			print ii
			print data[ii]
		ii=ii+1
	ii=0
	while ii < len(data):
		if (data[ii]['seq_end'] - data[ii]['seq_beg']+1) <= 5 and data[ii]['seq'] != 'seq': 
			print data
			data=np.delete(data,ii,0)
			print "###################################"
			print "PROBLEM: short monomer"
			print ii
			print data[ii]
		ii=ii+1
	data.sort(order='seq_beg')
	return data


################################
def re_estimation_monomer_limits(mono1, mono2, seq_file, mono_file):
	#############################
	#### temporary_file
	align_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	mono_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	if mono1['seq_beg'] < mono2['seq_end']:
		os.system("fastacmd -s "+mono2['seq']+" -L"+str(mono1['seq_beg'])+","+str(mono2['seq_end'])+" -d "+seq_file+">"+seq_tmpfile)
	else:
		print "error: bla bla bla ... " + str(mono1['seq_beg']) + " " + str(mono2['seq_end'])
		print mono1
		print mono2
		return(new_end, new_beg)
		quit()
	if mono1['strand'] == 'plus':
		os.system("fastacmd -s "+mono1['monomer']+" -d "+mono_file+">"+mono_tmpfile)
	else:
		os.system("fastacmd -s "+mono1['monomer']+" -S 2 -d "+mono_file+">"+mono_tmpfile)
	if mono2['strand'] == 'plus':
		os.system("fastacmd -s "+mono2['monomer']+" -d "+mono_file+" |grep -v \">\">>"+mono_tmpfile)
	else:
		os.system("fastacmd -s "+mono2['monomer']+" -S 2 -d "+mono_file+" |grep -v \">\">>"+mono_tmpfile)


	os.system("needle  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 0.5 -aformat3 fasta 2>1 > /dev/null")

	f=open(align_tmpfile,"r")
	try:
		alignment = AlignIO.read(f, "fasta")
	except:
		print "##################################"
		print "ERROR: alignement"
		print mono1
		print mono2
	f.close()
	os.remove(align_tmpfile)
	os.remove(seq_tmpfile)
	os.remove(mono_tmpfile)

	p1=0
	p2=0
	s1=0
	while p1 < mono1['len_mono']:
		if alignment[0].seq[p2] != '-':
			p1+=1
		if alignment[1].seq[p2] != '-':
			s1+=1
		p2+=1
	new_end=mono1['seq_beg']+s1-1
	while p1 <= mono1['len_mono']:
		if alignment[0].seq[p2] != '-':
			p1+=1
		if alignment[1].seq[p2] != '-':
			s1+=1
		p2+=1
	if alignment[1].seq[p2-1] == '-':
		s1+=1
	new_beg=mono1['seq_beg']+s1-1
	return(new_end, new_beg)

#############################
def get_sequence_monomers(data,seq):
	return data[data['seq']==seq]


#############################
def worker(work_queue, done_queue, data, seq_file, mono_file, blast_maximal_overlap, blast_minimal_length, waterman_minimal_score, waterman_minimal_length):
 
    print "starting: "+ current_process().name
    sys.stdout.flush()
    for ss in iter(work_queue.get, 'STOP'):
	print "#########################################"
	print current_process().name+" : "+ss
        sys.stdout.flush()
 	nb_selected_monomers=0

	data_seq=get_sequence_monomers(data,ss)
	data_seq=np.insert(data_seq,0,('seq', 'mono', 0.0, blast_minimal_length+1, -1, -1, -1, -1, 1000000, 0, 0, 'na'),0 )
	
	if len(data_seq[:]['len_seq']) == 1:
		print "ERROR: no blast hit for ... "+ss+ " is it normal???"
		continue
		
	seq_monomers_data=get_nonoverlapping_monomers(data_seq, blast_maximal_overlap, blast_minimal_length)
	
	print  "//////////////////////////"
	print  " data after selection for overlap ..."
	print  "//////////////////////////"
	seq_monomers_data=check_monomers(seq_monomers_data)
	print seq_monomers_data
	print  "//////////////////////////"

	if len(seq_monomers_data) > 1:
		seq_monomers_data.sort(order='seq_beg')

	if len(seq_monomers_data) == 1:
		seq_monomers_data=insert_newhit_by_smithwaterman_in_empty_data(seq_monomers_data, ss, 1, data_seq[0]['len_seq'], seq_file, mono_file,waterman_minimal_score, waterman_minimal_length)
		print  "//////////////////////////"
		print  " data after insertion ..."
		print  "//////////////////////////"
		seq_monomers_data=check_monomers(seq_monomers_data)
		print seq_monomers_data
		print  "//////////////////////////"
	
	if len(seq_monomers_data) > 1:
		search_new_initial_hits=1
		while search_new_initial_hits == 1:

			smd_len=len(seq_monomers_data)
			seq_monomers_data=insert_newhit_by_smithwaterman(seq_monomers_data, -1, 1, seq_file, mono_file,waterman_minimal_score, waterman_minimal_length)
			if len(seq_monomers_data) == smd_len:
				search_new_initial_hits=0
			if seq_monomers_data[0]['seq_beg'] < blast_minimal_length:
				search_new_initial_hits=0

		print  "//////////////////////////"
		seq_monomers_data=check_monomers(seq_monomers_data)
		print  " data after initial hits ..."
		print  "//////////////////////////"
		print seq_monomers_data
		print  "//////////////////////////"


		search_new_terminal_hits=1
		while search_new_terminal_hits == 1:
			seq_monomers_data=check_monomers(seq_monomers_data)

			last=len(seq_monomers_data)-1

			seq_monomers_data=insert_newhit_by_smithwaterman(seq_monomers_data, last, -1, seq_file, mono_file,waterman_minimal_score, waterman_minimal_length)
			if (len(seq_monomers_data)-1) == last:
				search_new_terminal_hits=0
			if (seq_monomers_data[len(seq_monomers_data)-1]['len_seq']-seq_monomers_data[len(seq_monomers_data)-1]['seq_end']) < blast_minimal_length:
				search_new_initial_hits=0

		print  "//////////////////////////"
		print  " data after terminal hits ..."
		print  "//////////////////////////"
		print seq_monomers_data
		print  "//////////////////////////"


	#### first pass for long hits: re-estimation
	if len(seq_monomers_data) > 1:
		seq_monomers_data.sort(order='seq_beg')




	ii=1
	while ii  < (len(seq_monomers_data)-1):
		print seq_monomers_data[ii]
		print seq_monomers_data[ii+1]
	
		dist_seq = seq_monomers_data[ii+1]['seq_beg']-seq_monomers_data[ii]['seq_end']
		dist_mono=(seq_monomers_data[ii]['len_mono'] -seq_monomers_data[ii]['mono_end'])+seq_monomers_data[ii+1]['mono_beg']
		if dist_seq > 1 and (math.fabs(dist_mono-dist_seq) < args.max_dist_hits or dist_seq < args.max_dist_hits):
			new_end, new_beg=re_estimation_monomer_limits(seq_monomers_data[ii], seq_monomers_data[ii+1], seq_file, mono_file)
			seq_monomers_data[ii]['seq_end']=new_end
			seq_monomers_data[ii+1]['seq_beg']=new_beg

			if seq_monomers_data[ii]['seq_end'] == seq_monomers_data[ii+1]['seq_end']:
				seq_monomers_data=np.delete(seq_monomers_data,ii+1 , 0)

		else:
			length=len(seq_monomers_data)
			seq_monomers_data=insert_newhit_by_smithwaterman(seq_monomers_data, ii, ii+1, seq_file, mono_file,waterman_minimal_score, waterman_minimal_length)
			if len(seq_monomers_data) > length:
				ii-=1
		seq_monomers_data=check_monomers(seq_monomers_data)
		seq_monomers_data.sort(order='seq_beg')
		ii+=1
	if len(seq_monomers_data) > 1:
		done_queue.put(seq_monomers_data)
	seq_monomers_data=check_monomers(seq_monomers_data)
    return True


##########################################
def writer_to_file(done_queue, seq_file, output_file):    
    while 1:
        m = done_queue.get()
	print "##########################################"
	if m != 'STOP':
		if (len(m)-1) == 2:
			#print len(m),1, m[1], m[1]['seq_beg']-1,m[1]['len_seq']-m[1]['seq_end']-1 
			get_sequence(m[1], seq_file, output_file)
	
		if (len(m)-1) > 2:
			#print len(m),1, m[1], m[1]['seq_beg']-1,m[2]['seq_beg']-m[1]['seq_end']-1 
			get_sequence(m[1], seq_file, output_file)

		ii=2
		while ii  < (len(m)-1):
			#print len(m),ii, m[ii], m[ii]['seq_beg']-m[ii-1]['seq_end']-1,m[ii+1]['seq_beg']-m[ii]['seq_end']-1 
			get_sequence(m[ii], seq_file, output_file)
			ii+=1
		ii=len(m)-1
		if ii > 2:
			#print len(m),ii, m[ii], m[ii]['seq_beg']-m[ii-1]['seq_end']-1,m[ii]['len_seq']-m[ii]['seq_end']-1 
			get_sequence(m[len(m)-1], seq_file, output_file)
		print "#####"
	else:
		break
    return True


#############################
parser = argparse.ArgumentParser(description='Search for tandemly repeated monomers ...')
parser.add_argument('-b', dest="blast_file")
parser.add_argument('-s', dest="seq_file")
parser.add_argument('-m', dest="mono_file")
parser.add_argument('-o', dest="output_file")
parser.add_argument('-proc',  dest="processors",       type=int,         default=4)
parser.add_argument('-verb',  dest="verbose",                   default=10)
parser.add_argument('-bmo',  dest="blast_maximal_overlap",     type=int, default=5)
parser.add_argument('-bmd',  dest="max_dist_hits",             type=int, default=30)
parser.add_argument('-bml',  dest="blast_minimal_length",    type=int, default=100)
parser.add_argument('-wml',  dest="waterman_minimal_length", type=int, default=100)
parser.add_argument('-wms',  dest="waterman_minimal_score",  type=float, default=200)
parser.add_argument('-wgo',  dest="waterman_extension_gap_penality", type=float, default=10)
parser.add_argument('-wge',  dest="waterman_opening_gap_penality",   type=float, default=0.5)
args = parser.parse_args()


#############################
# creating blast databases ...
os.system("formatdb -i "+args.seq_file+" -p F -o")
os.system("formatdb -i "+args.mono_file+" -p F -o")

if os.path.exists(args.output_file):
    os.remove(args.output_file)

#############################
# blasting ...
print 'blasting ...'

blastcmd='blastall -p blastn -F F -v 100000 -b 100000 -m 8 -a 4 '
blastcmd=blastcmd+' -i '+args.mono_file
blastcmd=blastcmd+' -d '+args.seq_file
blastcmd=blastcmd+' -o '+args.blast_file
blastcmd=blastcmd+' 1>2 > blast.log'

#os.system(blastcmd)


#############################
# geting sequences length ...
print 'getting sequence length ...'
length2_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
length1_tmpfile=tempfile.NamedTemporaryFile(delete=False).name

infoseqcmd1='infoseq '+args.seq_file+' -filter -only -name -length |sort -k1,1 > '+length1_tmpfile
infoseqcmd2='infoseq '+args.mono_file+' -filter -only -name -length |sort -k1,1 > '+length2_tmpfile

joincmd='sort -k1,1 '+args.blast_file
joincmd=joincmd+' |join - '+length2_tmpfile
joincmd=joincmd+' |sort -k2,2 '
joincmd=joincmd+' |join -1 2 -2 1 - '+length1_tmpfile
joincmd=joincmd+" |awk \' {$1=$1;ok=0} $9 > $10 {tmp=$9;$9=$10;$10=tmp;print $0, \"minus\";ok=1}  $9 < $10 && ok == 0 {print $0, \"plus\"} \'  >  "+args.blast_file+'4'

#os.system(infoseqcmd1)
#os.system(infoseqcmd2)
#os.system(joincmd)

#############################
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

