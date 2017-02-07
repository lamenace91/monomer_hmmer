#! /usr/bin/python3.4
# -*- coding: utf-8 -*-


# /usr/local/bin/python3.4





import get_monomers_lib as gm
import time
import pickle
import shutil as sh
from flask import Flask
from flask import render_template
from flask import send_from_directory
from flask import request
from flask import flash
from flask import redirect
from flask import url_for
import jinja2
import json







myport=5001

# Initialization of the Flask application
app = Flask(__name__)

app.config.from_object(__name__) 

# Path to working directory
working_dir="static/working_dir/"

#version number
version=0.1
verbose=5

# output filenames
reference_hmm_file = "reference.hmm"
genome_hmmer_file  = "genome.hmmer"
parsed_hmmer_file  = "genome.parsed.hmmer"



dump_data_file       = "data.dump"
dump_parameter_file  = "parameters.dump"
dump_other_data_file = "other_data.dump"

# maximum age of old files (days)
max_old_files=7

def myformat(value):
    return("{0:0.2f}".format(float(value)))
def FormatDecimal(value_list):
	if value_list != None:
		return(list(map(myformat, value_list)))
	else:
		return(None)
jinja2.filters.FILTERS['FormatDecimal'] = FormatDecimal






##########
# ROUTES #
##########


@app.route('/')
def nsa(): 
    return render_template('get_monomers.html', version=version)

			
#~ @app.route('/exec_previous',  methods=['GET', 'POST']) 
#~ def myexec_previous():   
	#~ idd = request.form['idd'].rstrip()
	#~ print("__"+idd+"__")	
	#~ fullname_dir=working_dir+'/'+idd
	#~ with open(fullname_dir+"/"+dump_data_file, 'rb') as inputfilea:
		#~ data= pickle.load(inputfilea)	
	#~ with open(fullname_dir+"/"+dump_parameters_file, 'rb') as inputfileb:
		#~ parameters= pickle.load(inputfileb)	
	#~ with open(working_dir+"/"+idd+"/"+dump_other_data_file, 'rb') as inputfilec:
		#~ other_data= pickle.load(inputfilec)	
	#~ if parameters['no_skel_table'] == 0:
		#~ return render_template('image_nsa_5.html',
				#~ no_skel_table=0,
				#~ working_dir=other_data['fullname_dir'],
				#~ idd=other_data['shortname_dir'], 
				#~ signal_plots=other_data['all_signal_plot_files'],
				#~ dd=data, 
				#~ image_dir=other_data['image_dir'],
				#~ parameters=parameters,
				#~ peak_nb=list(range(1, other_data['max_peak_nb']+1)),
				#~ ps=other_data['ps_all'],
				#~ datafile=other_data['datafile'],
				#~ zipdatafile=other_data['zipdatafile'])	
	#~ else:
		#~ return render_template('image_nsa_5.html',
				#~ no_skel_table=1,
				#~ working_dir=other_data['fullname_dir'],
				#~ idd=other_data['shortname_dir'], 
				#~ signal_plots=other_data['all_signal_plot_files'],
				#~ dd=data, 
				#~ image_dir=other_data['image_dir'],
				#~ parameters=parameters,
				#~ peak_nb=list(range(1, other_data['max_peak_nb']+1)),
				#~ ps=other_data['ps_all'],
				#~ datafile=other_data['datafile'],
				#~ zipdatafile=other_data['zipdatafile'])	
		#~ 



@app.route('/exec',  methods=['GET', 'POST']) 
def myexec():   
	
	####################################################################
	### removing old file
	gm.delete_old_files(working_dir, n_days=max_old_files)

	### temporary directory
	shortname_dir, fullname_dir=gm.create_temp_folder(directory=working_dir)

	### form parameters

	parameters=gm.get_parameters(request.form)
	genome_file=gm.get_genome_file(request.files, directory=fullname_dir)	
	reference_file=gm.get_reference_file(request.files, directory=fullname_dir)

	####################################################################
	if verbose > 0:
		print("Formating the hmmer db ...")
		t1=time.time()	
	monomer_hmm_file = gm.format_hmmdb(working_dir+"/"+shortname_dir+"/"+reference_file, db_name=working_dir+"/"+shortname_dir+"/"+reference_hmm_file)
	if verbose > 2:
		print("     hmm db in %s (%f s)" % (monomer_hmm_file, time.time()-t1))

	####################################################################
	if verbose > 0:
		print("Searching the sequences ...")
		t1=time.time()
	hmmer_ouput_file = gm.search_hmm(working_dir+"/"+shortname_dir+"/"+genome_file, 
										working_dir+"/"+shortname_dir+"/"+reference_hmm_file,
										output_file=working_dir+"/"+shortname_dir+"/"+genome_hmmer_file)
	if verbose > 2:
		print("     nhmmscan output in %s (%f s)" % (hmmer_ouput_file, time.time()-t1))
		
	####################################################################		
	if verbose > 0:
		print("Parsing the hmmer output ...")
		t1=time.time()
	hmm_hits = gm.parse_hmm_ouput(working_dir+"/"+shortname_dir+"/"+genome_hmmer_file)
	if verbose > 2:
		print("     number of hits: %d (%f s)" % (len(hmm_hits), time.time()-t1))

	####################################################################		
	if verbose > 0:
		print("Writing the hmmer output ...")
		t1=time.time()
	gm.write_all_blocks_all_seq(hmm_hits, working_dir+"/"+shortname_dir+"/"+parsed_hmmer_file)
	if verbose > 2:
		print("     done (%f s)" % (time.time()-t1))


	####################################################################		
	#### create archive files for downloads
	if verbose > 0:
		print("Zipping the output files ...")
		t1=time.time()
	sh.make_archive(working_dir+"/"+"archive_"+shortname_dir, 'zip', fullname_dir)
	sh.move(working_dir+"/"+"archive_"+shortname_dir+".zip", fullname_dir)
	if verbose > 2:
		print("     done (%f s)" % (time.time()-t1))

	 
	#####################
	#### create dump files for internal use
	## dump hits
	with open(fullname_dir+"/"+dump_data_file, 'wb') as output:
		pickle.dump(hmm_hits, output, pickle.HIGHEST_PROTOCOL)
		
	## dump parameters
	with open(fullname_dir+"/"+dump_parameter_file, 'wb') as outputb:
		pickle.dump(parameters, outputb, pickle.HIGHEST_PROTOCOL)

	## dump information
	with open(fullname_dir+"/"+dump_other_data_file, 'wb') as output:
		pickle.dump({'fullname_dir': fullname_dir,
					'shortname_dir': shortname_dir,
					'zipdatafile': "archive_"+shortname_dir+".zip"
				}, output, pickle.HIGHEST_PROTOCOL)
				 

	return render_template('get_monomers_out.html', 
			working_dir=fullname_dir,
			idd=shortname_dir, 
#			signal_plots=all_signal_plot_files,
#			dd=data, 
#			image_dir=image_dir,
			parameters=parameters,
#			peak_nb=list(range(1, max_peak_nb+1)),
#			ps=ps_all,
#			datafile=skeleton_interpol_data_file,
			zipdatafile="archive_"+shortname_dir+".zip")	
					
    
@app.route('/static/working_dir/<subdirectory>/<filename>', methods=['GET'])
def download(subdirectory, filename):
	print("######################################")
	print(subdirectory)
	print("######################################")
	return send_from_directory("static/working_dir/" + subdirectory  + "/", filename)
   
################################################################
@app.route('/hit_list/<idd>', methods=['POST', 'GET'])
def hit_list(idd, verbose=5):
	print(idd)
	verbose=5
	if verbose > 0:
		print("Reading hmmer files ...")
		t1 = time.time()
	table_data = gm.read_hit_list(working_dir+"/"+idd+"/"+parsed_hmmer_file)
	if verbose > 0:
		print("######################################")
		print(idd)
		print("######################################")
	return(json.dumps(table_data))

################################################################
@app.route('/summary_table/<idd>', methods=['POST', 'GET'])
def summary_table(idd, verbose=2):

	with open(working_dir+"/"+idd+"/"+dump_summary_file, 'rb') as inputfile:
		summary_table= pickle.load(inputfile)
	if verbose > 0:
		print("######################################")
		print(summary_table)
		print("######################################")
	return(json.dumps(summary_table))


@app.route('/data/<image>', methods=['GET'])
def data(image):
	return ''' Live long and prosper. !!!'''
    
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=myport, debug=True)
