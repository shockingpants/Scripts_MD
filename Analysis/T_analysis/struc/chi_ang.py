#!/opt/local/bin/python2.7
# chi_ang.py
# Extracts chi_ang from gromacs file
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
from mpl_toolkits.mplot3d import Axes3D
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
####        Script Startup
###############################################
#============================
#		Save Load
#============================
mode=[] # 0 --> run only 00 --> run or load  1 --> run and save  2 --> load
filename='dat0.db' 
mode,filename = SM.saveload(mode, filename) # Run this for checks.
if mode == 2:
	ds=SM.get_set(filename) # Loads saved data
start_time=time.time() # Start recording time
#============================
#		Plotting
#============================
MPL=MM.mpl([],'PCA_graph.pdf')
DPI=150

###############################################
####		     Functions
###############################################
##{{{
##}}}
###############################################
####		       Classes
###############################################
##{{{
class MyClass(object):
	def __init__(self):
		pass
	def __call__(self):
		pass
	def __str__(self):
		pass	
	def __getattr__(self):
		pass
	def __setattr__(self):
		pass
	def __delattr__(self):
		pass
	def __getitem__(self):
		pass
	def __setitem__(self):
		pass
	def __delitem__(self):
		pass
##}}}
###############################################
#### 		  Defining Parameters
###############################################
##{{{ Param
print 'Defining parameters...'
if mode==0 or mode==1:
	# Mode 0 --> run
	param={}
	param['step']=0.1 # Step
	param['start']=0 # Start
	param['end']=50 # End
elif mode==1:
	ds.add_param(param)
elif mode==2:
	param=ds.param

# To unload the parameter out of the dictionary
for i in iter(param):
	vars()[i]=param[i]
time=arange(param['start'],param['end'],param['step'])

##}}}
###############################################
####    		Main scripts
###############################################
##{{{
#=====================================
#	      	Extracting Chi
#=====================================
##{{{ Loading Chi Angles
print 'Opening and extracing from chi_ang file...     ('+tim()+')'
mmode=SM.msaveload(0,mode)
#-----------------Calculations--------------------
if mmode==0:
	folder='chi_angles'
	if os.path.isfile(folder+'/chi.log'):
		print 'Assumes that chi angles have already been processed by gromacs.\nProceeding...'
	else:
		os.environ['FOLD']=folder
		os.system('./chi_ang.sh')
		
	name=genfromtxt(folder+'/chi_index.txt', delimiter=" ", dtype=('float, |S20'))
	for ind,(res,i) in enumerate(name):
		print "Loading {0:s}...".format(i)
		if ind == 0:
			chi_ang=genfromtxt(folder+'/'+i,skip_header=12, dtype=(float, float))[:,1]
		else:
			chi_ang=np.vstack((chi_ang,genfromtxt(folder+'/'+i,skip_header=12, dtype=(float, float))[:,1]))
	chi_resid=name['f0'] # resid of residues with chi

	if mode==1:
		print 'Saving ... '
		ds.chi_ang=chi_ang
		ds.chi_resid=resid
elif mmode==2:
	print 'Loading ... '
	chi_ang=ds.chi_ang
	resid=ds.chi_resid
#-----------Secondary Calculations----------------
#a+b
#----------------Plotting-------------------------
print 'Plotting my results..    ('+tim()+')'
#fig1=plt.figure(1,figsize=(11.69, 8.27),dpi=DPI)
#ax11=fig1.add_subplot(111)
#lin1,=ax11.plot(a,b)
#MPL.show()
##}}}

##}}}



