#!/opt/local/bin/python2.7
#split_pdb.py
# Splits pdb into its different chain, aka ligand and receptors
#=============== MODULES ======================
##{{{
#-----------------  Suppress stdout --------------------
import sys
class NullDevice():
    def write(self, s):
        pass
ori = sys.stdout  # keep a reference to STDOUT
sys.stdout = NullDevice()  # redirect the real STDOUT
#--------------------------------------------------------
from numpy import *
import numpy as np
import scipy
import scipy.stats as ss
import os
import pprint
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
import MDAnalysis
from strtup import *
import math
import shelve
import shelve_manager as SM
import mpl_manager as MM
import re
import time
import atexit
import argparse
#========= R-python interface ===============
import rpy2.robjects as ro
R=ro.r
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
R.library('bio3d')
R.library('optparse')
R.library('cluster')
R.library('ncdf')
R.library('class')
R.library('stats4')
R.library('energy')

sys.stdout = ori  # turn STDOUT back on
##}}}
#==============================================

###############################################
####       Setting Up Parser
###############################################
##{{{
#=====================================
#        Handling Error
#=====================================
##{{{
class MyParser(argparse.ArgumentParser):
    def error(self, message):
		if re.search(r'\.py$',sys.argv[0]): # This mean from shell
			sys.stderr.write('error: %s\n' % message)
			self.print_help()
			sys.exit(2)
		else:
			raise Exception('error: %s\n' % message)
##}}}
#=====================================
#        Setting up
#=====================================
##{{{
parser=MyParser(description='Splitting pdb into separate chains.',usage='%(prog)s 1ycr.pdb -o rec.pdb lig.pdb')
parser.add_argument('input',type=argparse.FileType('r'))
parser.add_argument('-o','--output','--out',nargs='*',default=None,help='If not defined, creates ***_1.pdb, ***_2.pdb etc.')
parser.add_argument('-n',action='store_true',help='Checks for number of chains in input file without output.')
parser.add_argument('--ow',action='store_true',help='Auto overwrite or not(default).')
p=parser.parse_args()
f=p.input # Assign file
prot=[]
prot_info=[]
##}}}
##}}}
###########################################
####	  Saving Protein Data
###########################################
##{{{
#====== Extract
print 'Extracting lines from pdb...'
for line in f:
	if line.startswith('ATOM') or line.startswith('TER') or line.startswith('END'):
		prot.append(line)
	else:
		prot_info.append(line)

#====== Reporting number of chains
print 'Calculating number of chains...'
num=1
for ind,line in enumerate(prot):
	if line.startswith('TER') and prot[ind+1].startswith('ATOM'):
		num += 1

if p.n: # Calls to check number of chain
	print 'There are {0:0.0f} separate chains.'.format(num)
	if p.output==None: # If only -n, terminate prematurely
		sys.exit()

#===================================
#     		Saving
#===================================
##{{{
print 'Saving...'
#------- Split input Name
input_name=f.name.split('.') # Splits 1ycr.pdb into ['1ycr','pdb'] for example


#------- Create Filenames
if p.output==None:
	# create filenames eg.1ycr_1.pdb, 1ycr_2.pdb
	filenames=[''.join([input_name[0],'_',str(j+1),'.',input_name[1]]) for j in xrange(num)] 
	test=raw_input('Create a default list of proteins '+str(filenames)+'?(y/n)\n')
	if test=='y':
		pass
	elif test=='n':
		raise Exception('Please specify output files.')
	else:
		raise Exception('(Jon)Wrong option selected.')
else:
	if len(p.output)!=num:
		sys.stderr.write("ERROR: %d targets does not match %d chains in pdb\n" % (len(p.output),num))
		raise Exception("(Jon) %d targets does not match %d chains in pdb\n" % (len(p.output),num))
	elif len(p.output)==num:
		filenames=p.output


#------- Check filenames for existence
for i in filenames:
	if os.path.isfile(i):
		if p.ow:
			print i+' exists. Overwriting...'
		else:
			check=raw_input('File {0:s} exist. Overwrite? (y/n) \n'.format(i))
			if check=='y':
				print 'Overwriting...'
			elif check=='n':
				raise Exception('(Jon) File exist, please choose another name.')
			else:
				raise Exception('(Jon) Wrong option')


#------- Save files
filenames.reverse() # Facilitate pop
file_save=open(filenames.pop(),'w')
for ind,line in enumerate(prot):
	if line.startswith('ATOM'):
		file_save.write(line)
	elif line.startswith('TER') and prot[ind+1].startswith('ATOM'):
		file_save.write('END')
		file_save.close()
		file_save=open(filenames.pop(),'w') # Save into next file
	elif line.startswith('END'):
		file_save.write(line)
		file_save.close()

##}}}
##}}}

		
		
		

		



