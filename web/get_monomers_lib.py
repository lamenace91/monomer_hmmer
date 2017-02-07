import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import tempfile
from itertools import compress
import time
import shutil as sh
from werkzeug import secure_filename







##############################################################
### create temporary directory
###############################################################
def create_temp_folder(directory, verbose=0):
	fullname_dir=tempfile.mkdtemp(dir=directory)
	shortname_dir=os.path.basename(fullname_dir)
	if verbose >  5:
		print("#####")
		print("Directory full name: %s" % fullname_dir)
		print("Directory short name: %s" % shortname_dir)	
	return((shortname_dir, fullname_dir))


##############################################################
### get parameter from form
###############################################################

def get_parameters(request_form):
	parameters={}
	parameters['index'] = int(request_form['index'])
	return(parameters)


##############################################################
### get genome file
###############################################################
def get_genome_file(request_files, directory, verbose=0):
	### parsing file name from the form
	genome_file=request_files['genome_file']
	genome_filename = secure_filename(genome_file.filename)	
	if genome_filename != "":
		genome_file.save(os.path.join(directory, 'genome.fst'))
	else:
		sh.copyfile("static/genome_example2.fst", os.path.join(directory, "genome.fst"))
	return("genome.fst")


##############################################################
### get reference file
###############################################################
def get_reference_file(request_files, directory, verbose=0):
	### parsing file name from the form
	reference_file=request_files['reference_file']
	reference_filename = secure_filename(reference_file.filename)	
	if reference_filename != "":
		reference_file.save(os.path.join(directory, 'reference.fst'))
	else:
		sh.copyfile("static/reference_double_example.fst", os.path.join(directory, "reference.fst"))
	return("reference.fst")

##############################################################
### delete old files
###########################################################
def delete_old_files(path, n_days, verbose=9):	
	nb = 0	
	now = time.time()
	for dd in os.listdir(path):
		if os.stat(os.path.join(path, dd)).st_mtime > (now - n_days * 24 * 60 *60):
			if os.path.isdir(os.path.join(path, dd)):
				if verbose > 8:
					print("     Removing %s ..." % os.path.join(path, dd))
				sh.rmtree(os.path.join(path, dd))
				nb = nb + 1	
	return(nb)

##############################################################
### format hmmer db 
##############################################################
def format_hmmdb(seq_file, db_name=None, hmmpress=True):
	name = "monomer"  # ne pas modifier
	if db_name == None:
		db_name = tempfile.NamedTemporaryFile().name

	# create command line
	command = "hmmbuild -n " + name + " " + db_name + " " + seq_file
	# running the formatting
	result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)

	# hmmpress if needed 
	if hmmpress == True:
		command = "hmmpress " + db_name
		result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
	return(db_name)
	
##############################################################
### search a sequence file with a hmmer db 
##############################################################
def search_hmm(seq_file, monomer_hmm_file, output_file=None, cut_ga=True, options="", prog="nhmmer"):
	if output_file == None:
		output_file = tempfile.NamedTemporaryFile().name
	# create command line
	command = "nhmmer -o " + output_file
	if cut_ga == True:
		command = command + " --cut_ga "
	command = command + " " + options + " "
	command = command + " " + monomer_hmm_file + " " + seq_file
	# running the search
	result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=True)
	return(output_file)

##############################################################
### search a sequence file with a hmmer db 
##############################################################
def extract_region(regions, seq_file, output_file=None):
	if output_file == None:
		output_file = tempfile.NamedTemporaryFile().name

	fasta_sequences = SeqIO.parse(open(seq_file),'fasta')
	with open(output_file) as out_file:
		for fasta in fasta_sequences:
			if fasta.id in regions:
				for ii in range(0, len(regions[fasta.id])):
					print (regions[fasta.id][ii]['begin'])
					print (regions[fasta.id][ii]['end'])
					print (regions[fasta.id][ii]['strand'])
					
	return(output_file)




##############################################################
### parse a hmmer output
##############################################################
def parse_hmm_ouput(hmmer_ouput_file):
	output = {}
	fffile = open(hmmer_ouput_file,'r')
	hit = None
	for line in fffile:
		if line.startswith("Query:"):
			model_name, model_length = re.split("[ \n]+", line)[1:3]	

		if line.startswith("//"):  # to get the last hit
				if hit != None:
					output[target_name].append(hit)
				
		if line.startswith(">>"):
				if hit != None:
					output[target_name].append(hit)
			
				target_name = re.split("[ \n]+", line)[1]
				if target_name not in output:
					 output[target_name] = []
				fffile.readline()
				fffile.readline()
				line = fffile.readline()
				hit = get_hit(line)
				hit['monomer_seq']=""
				hit['target_seq']=""
				
		if 	re.match(r' *[^ ]+ +[0-9]+ .+ [0-9]+', line):
				if re.split("[ ]+", line)[1] == "monomer":
					monomer_seq_ali = re.split(" +", line)[3]
					hit['monomer_seq']=hit['monomer_seq'] + monomer_seq_ali
					#print("monomer: %s" % monomer_seq_ali)
				else:
					target_seq_ali  = re.split(" +", line)[3]
					hit['target_seq']=hit['target_seq'] + target_seq_ali
					#print("target: %s" % target_seq_ali)				
	print(len(output))
	
	for target_name in 	output.keys():
			output[target_name] = sorted(output[target_name], key=lambda k: k['alifrom'])		
			output[target_name] = remove_overlaping_hits_plusminus(output[target_name])
			output[target_name] = split_blocks(output[target_name], target_name)
			output[target_name] = reverse_all_minus_blocks(output[target_name])
			output[target_name] = junction_all_blocks(output[target_name])
			#output[target_name] = index_all_blocks(output[target_name])					
			
			#print_all_blocks(output[target_name])
	return(output)
	
##############################################################
### check and adjust monomer junctions
### all blocks
##############################################################
def junction_all_blocks(blocks):
	for ii in range(0, len(blocks)):
			blocks[ii]['index']=junction_one_blocks(blocks[ii])

	return(blocks)
##############################################################
### check and adjust monomer junctions
### one block
##############################################################
def junction_one_blocks(block):
		block['hits'][0]['junction']='first'
		for ii in range(1, block['nb_hits']):
			if block['hits'][ii]['envfromB'] == (block['hits'][ii-1]['envtoB'] + 1):
				block['hits'][ii]['junction'] = 'perfect_env-env'
			else:
				if block['hits'][ii]['envfromB'] == (block['hits'][ii-1]['alitoB'] + 1):
					block['hits'][ii]['junction'] = 'perfect_ali-env'
				else:
					if block['hits'][ii]['alifromB'] == (block['hits'][ii-1]['envtoB'] + 1):
						block['hits'][ii]['junction'] = 'perfect_env-ali'
					else:
						if block['hits'][ii]['alifromB'] == (block['hits'][ii-1]['alitoB'] + 1):
							block['hits'][ii]['junction'] = 'perfect_ali-ali'
						else:
							if block['hits'][ii]['envfromB'] > (block['hits'][ii-1]['envtoB'] + 1):
								block['hits'][ii]['junction'] = 'gap'
							else:
								block['hits'][ii]['junction'] = 'overlap'
		return(block)
##############################################################
### reverse all minus blocks
##############################################################
def reverse_all_minus_blocks(blocks):
	for ii in range(0, len(blocks)):
			blocks[ii]['index']=reverse_one_minus_block(blocks[ii])
	return(blocks)
	
##############################################################
### reverse one minus block
##############################################################
def reverse_one_minus_block(block):
	if block['hits'][0]['strand'] == 'plus':
		return(block)
	for ii in range(0, block['nb_hits']):
		print(block['length'])
		tmp                           = block['length'] - block['hits'][ii]['alitoB'] + 1
		block['hits'][ii]['alitoB']   = block['length'] - block['hits'][ii]['alifromB'] + 1
		block['hits'][ii]['alifromB'] = tmp
		tmp                           = block['length'] - block['hits'][ii]['envtoB'] + 1
		block['hits'][ii]['envtoB']   = block['length'] - block['hits'][ii]['envfromB'] + 1
		block['hits'][ii]['envfromB'] = tmp
		block['hits'][ii]['alifrom'] = block['end'] - block['hits'][ii]['alifrom'] + 1
		block['hits'][ii]['alito']   = block['end'] - block['hits'][ii]['alito'] + 1
		block['hits'][ii]['envfrom'] = block['end'] - block['hits'][ii]['envfrom'] + 1
		block['hits'][ii]['envto']   = block['end'] - block['hits'][ii]['envto'] + 1
	block['hits']=list(reversed(block['hits']))
	return(block)
	
##############################################################
### index all blocks
##############################################################
def index_all_blocks(blocks):
	for ii in range(0, len(blocks)):
			blocks[ii]['index']=index_one_block(blocks[ii])
	return(blocks)
	
##############################################################
### index one block
##############################################################
def index_one_block(block):
	output = [0] * (block['length']+1)
	for ii in range(1, block['hits'][0]['alifromB']):
		output[ii] = -1
	for kk in range(0, block['nb_hits']):
		ii=block['hits'][kk]['alifromB']
		index = block['hits'][kk]['hmmfrom']
		for jj in range (1, len(block['hits'][kk]['monomer_seq'])):
			if block['hits'][kk]['monomer_seq'][jj] == '.':
				output[ii] = -2
				ii = ii + 1
			else:
				if block['hits'][kk]['target_seq'][jj] != '-':
					output[ii] = index	
					ii = ii + 1
				#~ else:
					#~ print("%d %d - %c" %(jj, ii, 		block['hits'][kk]['target_seq'][jj]))
				index=index+1
				
		if kk < (block['nb_hits']-1):
			if (block['hits'][kk+1]['envfromB'] - block['hits'][kk]['envtoB']) == 1:
				# envfromB - envtoB link
				while ii <= block['hits'][kk]['envtoB']:
					output[ii] = index
					index = index + 1
					ii = ii + 1
				index = 1
				while ii < block['hits'][kk+1]['alifromB']:
					output[ii] = index
					index = index + 1
					ii = ii + 1
	for ii in range(block['hits'][block['nb_hits']-1]['alitoB'], block['length']):
		output[ii] = -1
	return(output)
##############################################################
### print all blocks
##############################################################
def print_all_blocks(blocks):
	for ii in range(0, len(blocks)):
		print(blocks[ii]['nb_hits'])
		for jj in 	range(0, blocks[ii]['nb_hits']):
			print ("block: %d / %d " % (ii+1,  len(blocks)))
			print ("block: %d - %d (%d bp)  " % (blocks[ii]['begin'],  blocks[ii]['end'], blocks[ii]['length']))
			print ("hit: %d / %d"   % (jj+1, blocks[ii]['nb_hits']))
			print ("ini ->   coord: %d - %d - %s" % (blocks[ii]['hits'][jj]['envfrom'],  blocks[ii]['hits'][jj]['envto'],  blocks[ii]['hits'][jj]["strand"]))
			print ("env ->   coord: %d - %d     " % (blocks[ii]['hits'][jj]['envfromB'], blocks[ii]['hits'][jj]['envtoB']))
			print ("ali ->   coord: %d - %d     " % (blocks[ii]['hits'][jj]['alifromB'], blocks[ii]['hits'][jj]['alitoB']))
			print ("%s" % (blocks[ii]['hits'][jj]['monomer_seq']))
			print ("%s" % (blocks[ii]['hits'][jj]['target_seq']))
		#~ alito=hit['alito']print_hit(blocks|ii][jj]
##############################################################
### print all blocks
##############################################################
def print_all_blocks3(blocks):
	for ii in range(0, len(blocks)):
		print(blocks[ii]['nb_hits'])
		for jj in 	range(0, blocks[ii]['nb_hits']):
			print ("%d / %d      %d - %d (%d bp)      %d / %d     %d - %d . %s %d - %d %d - %d  %s  %s" % (ii+1,  len(blocks), blocks[ii]['begin'],  blocks[ii]['end'], blocks[ii]['length'], jj+1, blocks[ii]['nb_hits'], blocks[ii]['hits'][jj]['envfrom'],  blocks[ii]['hits'][jj]['envto'],  blocks[ii]['hits'][jj]["strand"], blocks[ii]['hits'][jj]['envfromB'], blocks[ii]['hits'][jj]['envtoB'], blocks[ii]['hits'][jj]['alifromB'], blocks[ii]['hits'][jj]['alitoB'], blocks[ii]['hits'][jj]['monomer_seq'], blocks[ii]['hits'][jj]['target_seq']))
		#~ alito=hit['alito']print_hit(blocks|ii][jj]

##############################################################
### write all blocks into a file
##############################################################
def write_all_blocks(blocks, outfile, seqname=None):
	exist=False
	if  os.path.isfile(outfile): 
		exist=True
	handle=open(outfile,'a')
	if not exist:
		if seqname != None:
			handle.write("sequence ")
		handle.write("block block_nb block_begin block_end block_length")
		handle.write(" hit hit_nb envfrom envto strand envfromB envtoB alifromB alitoB junction")
		handle.write(" reference_seq monomer_seq\n")
	
	for ii in range(0, len(blocks)):
		for jj in 	range(0, blocks[ii]['nb_hits']):
			if seqname != None:
				handle.write(seqname+" ")			
			handle.write("%d %d %d %d %d %d %d %d %d %s %d %d %d %d %s %s %s\n" % (ii+1,  len(blocks), blocks[ii]['begin'],  blocks[ii]['end'], blocks[ii]['length'], jj+1, blocks[ii]['nb_hits'], blocks[ii]['hits'][jj]['envfrom'],  blocks[ii]['hits'][jj]['envto'],  blocks[ii]['hits'][jj]["strand"], blocks[ii]['hits'][jj]['envfromB'], blocks[ii]['hits'][jj]['envtoB'], blocks[ii]['hits'][jj]['alifromB'], blocks[ii]['hits'][jj]['alitoB'],  blocks[ii]['hits'][jj]['junction'], blocks[ii]['hits'][jj]['monomer_seq'], blocks[ii]['hits'][jj]['target_seq']))
	handle.close()
	return(0)	

##############################################################
### write all blocks for all sequences into a file
##############################################################
def write_all_blocks_all_seq(output, outfile):
	for target_name in 	output.keys():
		write_all_blocks(output[target_name], outfile, target_name)
	return(0)
	
##############################################################
### print all blocks
##############################################################
def print_all_blocks2(blocks):
	for ii in range(0, len(blocks)):
		print_one_block(blocks[ii])
		
##############################################################
### print one block
##############################################################
def print_one_block(block):
		print("########################")
		print("Nb of hits: %d" % block['nb_hits'])
		print("From-To: %d-%d" % (block['begin'],block['end']))
		print("Length: %d"     % block['length'])
		for jj in 	range(0, block['nb_hits']):
			print("#####")
			print ("hit: %d / %d      coord: %d - %d - %s" % (jj+1, block['nb_hits'], block['hits'][jj]['envfrom'], block['hits'][jj]['envto'], block['hits'][jj]["strand"]))
			print ("                  coord: %d - %d"      % (block['hits'][jj]['envfromB'], block['hits'][jj]['envtoB']))
			print ("                  coord: %d - %d"      % (block['hits'][jj]['alifromB'], block['hits'][jj]['alitoB']))
			print ("       %s" % (block['hits'][jj]['monomer_seq']))
			print ("       %s" % (block['hits'][jj]['target_seq']))
		#~ alito=hit['alito']print_hit(blocks|ii][jj]

##############################################################
### read hit list from tabular file
##############################################################
def read_hit_list(infile, sep=" "):
	out=[]
	with open(infile, "r") as f: 
		nb=0
		for line in f.readlines():
			li = line.lstrip()
			if not li.startswith("#"):
				nb=nb+1
				word = line.split(sep)
				if nb == 1:
					nbc=0
					name={}
					while nbc < len(word) and word[nbc] != '1':
						name[nbc]=word[nbc]
						nbc=nbc+1
				if nb > 1:
					tmp={}
					for ii in range(nbc):
						tmp[name[ii]]=word[ii]
					out.append(tmp)
	#print(out)
	return(out)


##############################################################
### get hit from hmmer output
##############################################################
def get_hit(line):				
	parsed_line = re.split("[ \n]+", line)
	nb=len(parsed_line)
	hit = {}
	hit["score"] = float(parsed_line[2])
	hit["evalue"] = float(parsed_line[4])
	hit["strand"] = "plus"
	hit["hmmfrom"] = int(parsed_line[5])
	hit["hmmto"] = int(parsed_line[6])
	hit["strand"] = "plus"
	hit["envfrom"] = int(parsed_line[11])
	hit["envto"] = int(parsed_line[12])
	hit["alifrom"] = int(parsed_line[8])
	hit["alito"] = int(parsed_line[9])
	if hit["alifrom"] > hit["alito"]:
		hit["strand"] = "minus"
		hit["alifrom"] = int(parsed_line[9])
		hit["alito"] = int(parsed_line[8])
		hit["envfrom"] = int(parsed_line[12])
		hit["envto"] = int(parsed_line[11])
	return(hit)					
##############################################################
### get hit from hmmer output
##############################################################
def get_hit2(file_iterator):
	file_iterator.newline()
	file_iterator.newline()
	line = file_iterator.newline()
	parsed_line = re.split("[ \n]+", line)
	nb=len(parsed_line)
	hit = {}
	hit["score"] = float(parsed_line[2])
	hit["evalue"] = float(parsed_line[4])
	hit["strand"] = "plus"
	hit["hmmfrom"] = int(parsed_line[5])
	hit["hmmto"] = int(parsed_line[6])
	hit["strand"] = "plus"
	hit["envfrom"] = int(parsed_line[11])
	hit["envto"] = int(parsed_line[12])
	hit["alifrom"] = int(parsed_line[8])
	hit["alito"] = int(parsed_line[9])
	if hit["alifrom"] > hit["alito"]:
		hit["strand"] = "minus"
		hit["alifrom"] = int(parsed_line[9])
		hit["alito"] = int(parsed_line[8])
		hit["envfrom"] = int(parsed_line[12])
		hit["envto"] = int(parsed_line[11])
	return(hit)					
	
		
##############################################################
### remove overlaping hits on different strands
##############################################################
def remove_overlaping_hits_plusminus(hits):
	nb_hits = len(hits)
	as_to_be_conserved = [True] * nb_hits
	for ii in range(0, nb_hits): 
		as_to_be_conserved[ii]= is_best_overlapping_hit(hits, ii)
	#print(as_to_be_conserved)
	output = list(compress(hits, as_to_be_conserved))
	return(output)
##############################################################
### compare and return the best hit
##############################################################
def is_best_overlapping_hit(hits, hit_index, overlap=20):
	nb_hits = len(hits)
	ii = hit_index+1
	output = True
	while output == True and ii < nb_hits  and (hits[ii]['envfrom']-overlap)  <= hits[hit_index]['envto']:
		if hits[ii]['strand'] != hits[hit_index]['strand'] and hits[ii]['score'] > hits[hit_index]['score']:
			output = False
		ii = ii+1
	return(output)
	
##############################################################
### split hmmer output into blocks
##############################################################
def split_blocks(hits, seqname,distance=50):
	nb_hits = len(hits)
	output = []

	block = {}
	block['sequence' ]= seqname
	block['hits']     = []
	block['hits'].append(hits[0])
	block['begin']   = hits[0]['envfrom']
	block['end']     = hits[0]['envto']
	block['nb_hits'] = 1
#	print(block)
	for ii in range(1, nb_hits): 
		# not on the same strand
		# or
		# separated by more than distance
		if hits[ii-1]['strand'] != hits[ii]['strand']  or  hits[ii-1]['envto'] < (hits[ii]['envfrom'] - distance):
			block['length'] = block['end'] - block['begin'] + 1
			block = get_hit_positions_along_block(block)
			output.append(block)
			block = {}
			block['hits'] = []
			block['hits'].append(hits[ii])
			block['begin'] = hits[ii]['envfrom']
			block['end'] = hits[ii]['envto']
			block['nb_hits'] = 1
		else:
			block['hits'].append(hits[ii])
			block['end'] = hits[ii]['envto']
			block['nb_hits'] = block['nb_hits'] + 1
	block['length'] = block['end'] - block['begin'] + 1
	block = get_hit_positions_along_block(block)
	output.append(block)

	return(output)
	
##############################################################
### calculate hit positions along the blocks instead of 
###     positions along the sequence
##############################################################
def get_hit_positions_along_block(block):
	output=block
	for ii in range(0, block['nb_hits']):
		block['hits'][ii]['alifromB'] = block['hits'][ii]['alifrom'] - block['begin'] + 1
		block['hits'][ii]['alitoB']   = block['hits'][ii]['alito']   - block['begin'] + 1
		block['hits'][ii]['envfromB'] = block['hits'][ii]['envfrom'] - block['begin'] + 1
		block['hits'][ii]['envtoB']   = block['hits'][ii]['envto']   - block['begin'] + 1	
	return(block)

#
