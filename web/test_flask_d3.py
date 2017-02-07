#! /usr/bin/python3.4
# -*- coding: utf-8 -*-


# /usr/local/bin/python3.4




import random
from flask import Flask
from flask import render_template
from flask import send_from_directory
from flask import request
from flask import flash
from flask import redirect
from flask import url_for
import jinja2
import json
import pandas
import numpy as np






myport=5001

# Initialization of the Flask application
app = Flask(__name__)

app.config.from_object(__name__) 






##########
# ROUTES #
##########

  
################################################################
@app.route('/<nbseqmax>/<nbasmax>', methods=['POST', 'GET'])
def generate_data(nbseqmax, nbasmax):
	nbseq=random.randint(1, int(nbseqmax) )
	maxas=0
	maxposition=0
	data=[]
	#~ begin=[]
	#~ end=[]
	#~ seq=[]
	for ii in range(0,nbseq):
		dd={"name": "seq_"+str(ii), "alphasat": []}
		position=random.randint(1, 200 )
		nbas=random.randint(1, int(nbasmax))
		if nbas > maxas:
			maxas=nbas
		for jj in range(0,nbas):
			length=random.randint(100, 200 )
			dd["alphasat"].append({"begin": position, "end": position+length, "seq": ii, "strand": random.randint(0,1)})
			position=position+length+random.randint(1, 200 )
			if position > maxposition:
				maxposition=position
		data.append(dd)
	print(data)
	return render_template('test_flask_d3.html', data=data, nbseq=nbseq, maxas=maxas, maxposition=maxposition)
  
if __name__ == "__main__":
    app.run(host='0.0.0.0', port=myport, debug=True)
