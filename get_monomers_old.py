#! /usr/bin/python

import os
import csv
import numpy as np
from numpy import genfromtxt
from Bio import AlignIO
import time
import argparse
import tempfile
from Bio import SeqIO

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
	print "---------------------------------------"	
	print seq1
	print seq2
	print beg
	print end
	return (beg, end)
#############################################################
def get_similarity(seq1, seq2):
	print "///////////////////"
	simil=0.0
	total=0.0
	ii=0
	while ii < len(seq1):
		if seq1[ii] != '-' and seq2[ii] != '-':	
			if seq1[ii] == seq2[ii]:	
				simil+=1.0
			total+=1.0
		ii+=1
	print "///////////////////"
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
	os.system("fastacmd -s "+mono['seq']+" -L"+str(mono['seq_beg'])+","+str(mono['seq_end'])+" -d "+seq_file+">"+seq_tmpfile)
	for mono_record in SeqIO.parse(mono_file, "fasta"):
   	 	print mono_record.id
		if mono['strand'] == 'plus':
			os.system("fastacmd -s "+mono_record.id+" -d "+mono_file+">"+mono_tmpfile)
		else:
			os.system("fastacmd -s "+mono_record.id+" -S 2 -d "+mono_file+">"+mono_tmpfile)
		os.system("needle  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 0.5 -aformat3 fasta 2>1 > /dev/null")
		alignment = AlignIO.read(align_tmpfile, "fasta")		
		simil=get_similarity(alignment[0].seq,alignment[1].seq)
		if simil > simil_max:
			beg,end=get_limits_on_mono(alignment[0].seq,alignment[1].seq)
			new_mono['simil']=simil
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
#		print beg
#		print end
		if len(sub_data_tmp4) == 0:
			return new_mono	
#		print set(sub_data_tmp4['monomer'])
		
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
#		print beg
#		print end
		if len(sub_data_tmp4) == 0:
			return new_mono	
#		print set(sub_data_tmp4['monomer'])
		for mm in set(sub_data_tmp4['monomer']):
			sub_data_tmp5=sub_data_tmp4[sub_data_tmp4[:]['monomer']==mm]
			for ss in set(sub_data_tmp5['strand']):
				sub_data_tmp=sub_data_tmp5[sub_data_tmp5[:]['strand']==ss]
				ind = np.argmax(sub_data_tmp[:]['score'])
				new_mono=sub_data_tmp[ind]
				sub_data_tmp[ind]['score']=-1
				selected_monomers=1
				new_mono['score']=1
			
				print "###########"
				print mm
			
#				print "~~~~~~~~~~~~~~~~~~~~~~~"
#				print sub_data_tmp
#				print "~~~~~~~~~~~~~~~~~~~~~~~"
			
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
				#print "ssssssssssssssssssss"
				#print new_len
				#print new_mono['seq_end']-new_mono['seq_beg']
				#print new_mono2
#				print "---->   "+str(new_mono)
		
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
#		print sub_data_tmp[ind]
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
def search_nonoverlaped_completed_monomers(data,overlap_len):
	nb_selected_monomers = 0
	selected_monomers_data=[]
	selected_monomers=[]
	if len(data) == 0:
		return selected_monomers_data
	ind = np.argmax(data[:]['score'])
	while len(selected_monomers) < len(data) and data[ind]['score'] > 0:
		overlap=0
		for sm in  selected_monomers:
			if data[ind]['seq_beg'] < (data[sm]['seq_end']-overlap_len) and data[ind]['seq_end'] >  (data[sm]['seq_beg']+overlap_len):
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
		ind = np.argmax(data[:]['score'])
	if len(selected_monomers_data) > 1:
		selected_monomers_data.sort(order='seq_beg')
	return(selected_monomers_data)
#	fSubset = open("monomer_"+ss+".txt", 'w') 
#	np.savetxt(fSubset, selected_monomers_data,fmt="%s",delimiter=":")
#	fSubset.close()

################################
def re_estimation_monomer_limits(mono1, mono2, seq_file, mono_file):
	#############################
	#### temporary_file
	align_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	seq_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	mono_tmpfile=tempfile.NamedTemporaryFile(delete=False).name
	#align_tmpfile="needle_tmp"
	#seq_tmpfile="seq_tmp"
	#mono_tmpfile="mono_tmp"	
	#############################
	#### alignment
#	if   len(mono1) == 0:
#		# initial monomer
#		os.system("fastacmd -s "+mono2['seq']+" -L"+str(1)+","+str(mono2['seq_end'])+" -d "+seq_file+">"+seq_tmpfile)
#	elif len(mono2) == 0:
#		# terminal monomer
#		os.system("fastacmd -s "+mono1['seq']+" -L"+str(mono1['seq_beg'])+","+str(mono1['len_seq'])+" -d "+seq_file+">"+seq_tmpfile)
#	else:
#		# internal monomer
	os.system("fastacmd -s "+mono2['seq']+" -L"+str(mono1['seq_beg'])+","+str(mono2['seq_end'])+" -d "+seq_file+">"+seq_tmpfile)
	
	
#	if len(mono1) != 0:	
	if mono1['strand'] == 'plus':
		os.system("fastacmd -s "+mono1['monomer']+" -d "+mono_file+">"+mono_tmpfile)
	else:
		os.system("fastacmd -s "+mono1['monomer']+" -S 2 -d "+mono_file+">"+mono_tmpfile)


#	if len(mono2) != 0:	
	if mono2['strand'] == 'plus':
		os.system("fastacmd -s "+mono2['monomer']+" -d "+mono_file+" |grep -v \">\">>"+mono_tmpfile)
	else:
		os.system("fastacmd -s "+mono2['monomer']+" -S 2 -d "+mono_file+" |grep -v \">\">>"+mono_tmpfile)


	os.system("needle  -asequence "+mono_tmpfile+" -bsequence "+seq_tmpfile+" -outfile "+align_tmpfile+" -gapopen 10 -gapextend 0.5 -aformat3 fasta 2>1 > /dev/null")
	alignment = AlignIO.read(align_tmpfile, "fasta")		
	os.remove(align_tmpfile)
	os.remove(seq_tmpfile)
	os.remove(mono_tmpfile)
	
	
	#############################
	#### parsing
	#	if    len(mono1) == 0:
	#		######################################
	#		print "" 
	#		
		
	#	elif  len(mono2) == 0: 
	#		######################################
	#		print "" 
		
	#		
	#	else:
	#		######################################
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
#	print mono1['seq_beg']
#	print p1
#	print p2
#	print s1
	while p1 <= mono1['len_mono']:
		if alignment[0].seq[p2] != '-':	
			p1+=1
		if alignment[1].seq[p2] != '-':	
			s1+=1	
		p2+=1
	if alignment[1].seq[p2-1] == '-':
		s1+=1	
	new_beg=mono1['seq_beg']+s1-1
		
		
		
#	print alignment[0].seq
#	print alignment[1].seq
#	print alignment[1].seq[p2]
#	print p1
#	print p2
#	print s1
	
	return(new_end, new_beg)

#############################
parser = argparse.ArgumentParser(description='Search for tandemly repeated monomers ...')
parser.add_argument('-b', dest="blast_file")
parser.add_argument('-s', dest="seq_file")
parser.add_argument('-m', dest="mono_file")
parser.add_argument('-O', dest="output_file")

parser.add_argument('-o',  dest="overlap",        default=5)
parser.add_argument('-l',  dest="min_hit_len", default=100)
parser.add_argument('-x',  dest="max_dist_hits",  default=30)
parser.add_argument('-v',  dest="verbose",  default=10)
args = parser.parse_args()
#parser.print_help()

#blast_file='chrAll.ALR_Alpha.monomer.blastn2',
#seq_file="chrAll.ALR_Alpha.fst"
#mono_file="monomer_types.fst"
#min_hit_len=100
#overlap_length=5 # max overlap authorized for two blast hits
#limit_dist=10    # max distance authorized between sequence gap and monomer gap



#############################
# creating blast databases ...
os.system("formatdb -i "+args.seq_file+" -p F -o")
os.system("formatdb -i "+args.mono_file+" -p F -o")

if os.path.exists(args.output_file):
    os.remove(args.output_file)







#############################
# blasting ...
print 'blasting ...'

blastcmd='blastall -p blastn -F F -v 10000 -b 10000 -m 8 -a 4 '
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
joincmd=joincmd+" |awk \' {$1=$1} $9 > $10 {tmp=$9;$9=$10;$10=tmp; print $0, \"minus\"}  $9 < $10 {print $0, \"plus\"} \' > "+args.blast_file+'4'

#os.system(infoseqcmd1)
#os.system(infoseqcmd2)
#os.system(joincmd)

#############################
# reading blast output ...
data = np.genfromtxt(args.blast_file+'4', delimiter=' ', usecols=[0,1,2,3,6,7,8,9,11,12,13,14], 
      names=['seq','monomer','simil', 'length','mono_beg', 'mono_end', 'seq_beg', 'seq_end', 'score', 'len_mono', 'len_seq', 'strand']  , dtype=None)


print "begining ..."
subi=0
subj=subi
for ss in set(data['seq']):
	nb_selected_monomers=0
	sub_data1=data[data[:]['seq_beg'] < data[:]['seq_end']]
	sub_data2=sub_data1[sub_data1[:]['seq']==ss]
	sub_data=sub_data2[sub_data2[:]['length'] >=args.min_hit_len]

	selected_monomers_data=search_nonoverlaped_completed_monomers(sub_data, args.overlap)

#	print selected_monomers_data
	if len(selected_monomers_data) >= 1:
#		print selected_monomers_data[0]['seq_beg']
#		news_monomer=build_monomers_from_small_hits2(data, ss, 1,selected_monomers_data[0]['seq_beg'] , 5)
		begin_length=selected_monomers_data[0]['seq_beg']
		print "begin: "+str(begin_length)
		while begin_length > args.min_hit_len and begin_length > 0:
			news_monomer=get_monomers_from_small_hits(data, ss, 1,begin_length)
			if len(news_monomer) > 0:
				selected_monomers_data=np.append(selected_monomers_data,news_monomer)
#				print news_monomer
				print "   initial: insertion of new monomers: "+str(len(news_monomer))
				begin_length=news_monomer['seq_beg']
				print news_monomer
			else:
				print "   initial: insertion of new monomer failed"
				begin_length=-1
			selected_monomers_data.sort(order='seq_beg')
			print begin_length
			
	
#			news_monomer=build_monomers_from_small_hits2(data, ss, selected_monomers_data[len(selected_monomers_data)-1]['seq_end'],selected_monomers_data[len(selected_monomers_data)-1]['len_seq'] , 5)
		news_monomer=get_monomers_from_small_hits(data, ss, selected_monomers_data[len(selected_monomers_data)-1]['seq_end'],selected_monomers_data[len(selected_monomers_data)-1]['len_seq'])	
		if len(news_monomer) > 0:
#			print news_monomer
			selected_monomers_data=np.append(selected_monomers_data,news_monomer)
			print "   terminal: insertion of new monomers: "+str(len(news_monomer))
		else:
			print "   terminal: insertion of new monomer failed"
		selected_monomers_data.sort(order='seq_beg')
		
	if len(selected_monomers_data) > 1:
		############################
		## initial monomer
		#dist_seq=selected_monomers_data[0]['seq_beg']
		#dist_mono=selected_monomers_data[0]['mono_beg']
		#if dist_seq < args.min_hit_len and (((dist_mono-dist_seq) < args.max_dist_hits  and (dist_seq-dist_mono) < args.max_dist_hits) or dist_seq < args.max_dist_hits):
			#new_end, new_beg=re_estimation_monomer_limits([], selected_monomers_data[0], args.seq_file, args.mono_file)
			#selected_monomers_data[0]['seq_beg']=new_beg
	
		ii=0
		while ii  < (len(selected_monomers_data)-1):
			print "#########"
			print "ii:", str(ii)
			print "len: "+str(len(selected_monomers_data))
			if args.verbose  > 2:
				print selected_monomers_data
				
			if args.verbose  > 0:
				print selected_monomers_data[ii]
				print selected_monomers_data[ii+1]
		
			dist_seq = selected_monomers_data[ii+1]['seq_beg']-selected_monomers_data[ii]['seq_end']
			dist_mono=(selected_monomers_data[ii]['len_mono'] -selected_monomers_data[ii]['mono_end'])+selected_monomers_data[ii+1]['mono_beg']
			
			if dist_seq == 1: 
				print "   junction: perfect"
			else:
				if dist_seq > 1 or dist_seq <= 0:
					#if dist_seq < args.min_hit_len and (((dist_mono-dist_seq) < args.max_dist_hits  and (dist_seq-dist_mono) < args.max_dist_hits) or dist_seq < args.max_dist_hits):
					print "#############"
					print dist_mono
					print dist_seq
					print args.max_dist_hits
					if (((dist_mono-dist_seq) < args.max_dist_hits  and (dist_seq-dist_mono) < args.max_dist_hits) or dist_seq < args.max_dist_hits):
						#print ("alignment...")
						# to extract monomers
						#to extract target sequence
						#print selected_monomers_data[ii]['seq_beg']
						#print selected_monomers_data[ii+1]['seq_end']
						new_end, new_beg=re_estimation_monomer_limits(selected_monomers_data[ii], selected_monomers_data[ii+1], args.seq_file, args.mono_file)
				
						selected_monomers_data[ii]['seq_end']=new_end
						selected_monomers_data[ii+1]['seq_beg']=new_beg
						if selected_monomers_data[ii]['seq_end'] == selected_monomers_data[ii+1]['seq_end']:
							selected_monomers_data=np.delete(selected_monomers_data,ii+1 , 0)
						if new_beg-new_end == 1:
						
							print "   junction: restimation perfect"
						else:
							print ("   junction: restimation imperfect"+str(new_beg-new_end))
						print selected_monomers_data[ii]
						print selected_monomers_data[ii+1]
					else:
						
						region_beg=selected_monomers_data[ii]['seq_end']
						region_end=selected_monomers_data[ii+1]['seq_beg']
			#			news_monomer=build_monomers_from_small_hits2(data, ss, region_beg, region_end, 5)
						news_monomer=get_monomers_from_small_hits(data, ss, region_beg, region_end)
						print "-------------> "+str(news_monomer)
						if len(news_monomer) > 0:
							selected_monomers_data=np.append(selected_monomers_data,news_monomer)
							ii=ii-1
							print "   junction: insertion of new monomers: "+str(len(news_monomer))
						else:
							print "   junction: insertion of new monomer: failed"
						#print len(selected_monomers_data)
						selected_monomers_data.sort(order='seq_beg')
						
			ii+=1				
			

		############################
		## terminal monomer
#		last=len(selected_monomers_data)-1
#		dist_seq=selected_monomers_data[last]['len_seq']-selected_monomers_data[last]['seq_end']
#		dist_mono=selected_monomers_data[last]['len_mono']-selected_monomers_data[last]['mono_end']
#		if dist_seq < args.min_hit_len and (((dist_mono-dist_seq) < args.max_dist_hits  and (dist_seq-dist_mono) < args.max_dist_hits) or dist_seq < args.max_dist_hits):
#			new_end=re_estimation_terminal_monomer_limit(selected_monomers_data[last], args.seq_file, args.mono_file)
#			selected_monomers_data[last]['seq_end']=new_end
		

		ii=0
		while ii  < (len(selected_monomers_data)):
			mono=selected_monomers_data[ii]['seq_end']-selected_monomers_data[ii]['seq_beg']+1
			if mono < 70:
				selected_monomers_data=np.delete(selected_monomers_data, ii, 0)
			ii+=1
			
		ii=0
		while ii  < (len(selected_monomers_data)-1):
			dist_seq = selected_monomers_data[ii+1]['seq_beg']-selected_monomers_data[ii]['seq_end']
			if dist_seq != 1 and dist_seq <  args.min_hit_len:
				new_end, new_beg=re_estimation_monomer_limits(selected_monomers_data[ii], selected_monomers_data[ii+1], args.seq_file, args.mono_file)
				selected_monomers_data[ii]['seq_end']=new_end
				selected_monomers_data[ii+1]['seq_beg']=new_beg
				if selected_monomers_data[ii]['seq_end'] == selected_monomers_data[ii+1]['seq_end']:
					selected_monomers_data=np.delete(selected_monomers_data,ii+1 , 0)

			
			ii+=1
	
		
		ii=0
		while ii  <= (len(selected_monomers_data)-1):
			mono=selected_monomers_data[ii]['seq_end']-selected_monomers_data[ii]['seq_beg']+1
			
			if selected_monomers_data[ii]['seq_beg'] < selected_monomers_data[ii]['seq_end'] and mono > 10:
				new_mono=align_mono(selected_monomers_data[ii], args.seq_file, args.mono_file)
				selected_monomers_data[ii]=new_mono
				if ii == 0:
					before=selected_monomers_data[ii]['seq_beg']-1
					after=selected_monomers_data[ii+1]['seq_beg']-selected_monomers_data[ii]['seq_end']-1
				elif ii==(len(selected_monomers_data)-1):
					before=selected_monomers_data[ii]['seq_beg']-selected_monomers_data[ii-1]['seq_end']-1
					after=selected_monomers_data[ii]['len_seq']-selected_monomers_data[ii]['seq_end']
				else:	
					before=selected_monomers_data[ii]['seq_beg']-selected_monomers_data[ii-1]['seq_end']-1
					after=selected_monomers_data[ii+1]['seq_beg']-selected_monomers_data[ii]['seq_end']-1

				print (str(ii)+" _monomer_: "+str(before)+" "+str(mono)+" "+str(after)+" "+str(selected_monomers_data[ii]))
				os.system("fastacmd -s "+selected_monomers_data[ii]['seq']+" -L "+str(selected_monomers_data[ii]['seq_beg'])+","+str(selected_monomers_data[ii]['seq_end'])+" -S 1 -d "+args.seq_file+" >> "+ args.output_file)
			else:
				selected_monomers_data=np.delete(selected_monomers_data,ii, 0)
				ii=ii-1
			ii+=1				
