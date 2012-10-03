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
filename='chi_ang0.db' 
mode,filename = SM.saveload(mode, filename) # Run this for checks.
if mode == 2 or mode == 1 or mode==0.1:
	ds=SM.get_set(filename) # Loads saved data
start_time=time.time() # Start recording time
#============================
#		Plotting
#============================
MPL=MM.mpl([],'CCA_graph.pdf')
DPI=150

###############################################
####		     Functions
###############################################
##{{{
def distance(coord,ref,n=3,t=[]):
	"""Calculates the distance using the xyz mat and a reference structure (slow)"""
	import time
	s_time=time.time()
	if np.shape(coord)[0] != len(ref):
		# Checks Dimension between coord and ref
		print np.shape(coord)
		print len(ref)
		raise Exception('Dimensions do not match. Make sure its n*p, 1*n')
	if np.shape(coord)[0]%n!=0:
		# Check dimension against n
		print np.shape(coord)[0]
		raise Exception('This only works with system with {0:0.0f} dimension.'.format(n))
	temp,ave_mat=np.meshgrid(np.arange(np.shape(coord)[1]),ref)
	diff_mat=coord-ave_mat # Difference matrix
	diff_mat2=diff_mat**2 # Difference matrix squared
	sum_mat=[]
	for i in arange(len(ref))[::n]:
		sum_=np.sum(diff_mat2[i:i+n,:],axis=0)
		if i==0:
			sum_mat=sum_
		else:
			sum_mat=vstack((sum_mat,sum_))
	if t != []:
		print time.time()-s_time
	return np.sqrt(sum_mat)



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
if mode==0 or mode==1 :
	# Mode 0 --> run
	param={}
	param['step']=0.1 # Step
	param['start']=0 # Start
	param['end']=50 # End
	if mode==1:
		ds.add_param(param)
elif mode==2 or mode==0.1:
	ds.get_param()
	
##}}}
###############################################
####    		Main scripts
###############################################
##{{{
#=====================================
#	      	Extracting Chi
#=====================================
##{{{ Loading Chi Angles
print 'Opening and extracting from chi_ang file...     ('+tim()+')'
with SM.sl(['chi_ang','chi_resid'],ds,2,mode):
	folder='chi_angles' # Folder that stores the temporary chi files
	# Process shell script
	if os.path.isfile(folder+'/chi.log') and os.path.isfile(folder+'/chi_index.txt'):
		print 'chi_index.txt exist. Assumes that chi angles have already been processed by gromacs.\nProceeding...'
	else:
		os.environ['FOLD']=folder
		os.system('./chi_ang.sh')

	# Extract data from chi files
	name=genfromtxt(folder+'/chi_index.txt', delimiter=" ", dtype=('|S20, |S20'))
	for ind,(res,i) in enumerate(name):
		print "Loading {0:s}...".format(i)
		if ind == 0:
			chi_ang=genfromtxt(folder+'/'+i,skip_header=12, dtype=(float, float))[:,1]
		else:
			chi_ang=np.vstack((chi_ang,genfromtxt(folder+'/'+i,skip_header=12, dtype=(float, float))[:,1]))
	chi_resid=name['f0'] # resid of residues with chi
	
	# Remove files
	import shutil
	shutil.rmtree(folder) #remove Chi Angles
	#remove extra trajectory
	try:
		os.remove('../../prod_vac.trr')
	except OSError:
		filen=raw_input('prod_vac.trr not found. What is the name of the temporary gromacs trajectory?\n')
		try:
			os.remove('../../'+filen)
		except OSError:
			print 'Unable to remove temporary file, please remove it manually.'

print '{0:0.1f}s has elapsed.'.format(time.time()-start_time)
##}}}
#=====================================
#	      	Extracting xyz
#=====================================
##{{{
print 'Loading from structural file...     ('+tim()+')'
#---------------------------------------------------------
struc_file='dat0.db'
struc=SM.get_set(struc_file) #import data shelve
b3d=importr('bio3d') #import bio3d
with SM.sl(['xyz','fluc','ave','cij'],struc,2,2):
	pass
# Calculating distance from avg struture
xyz2=np.array(xyz)
xyz2=xyz2.T
with SM.sl(['D'],ds,2,mode):
	D=distance(xyz2,np.array(ave),t='y') #D of RMDS. Distance from ave structure
print '{0:0.1f}s has elapsed.'.format(time.time()-start_time)
##}}}
#=====================================
#	Canonical Correlation Analysis
#=====================================
##{{{ Performing CCA
print 'Performing CCA...     ('+tim()+')'
#-----------------Calculations--------------------
#with SM.sl([],ds,1,mode):
CCA=importr('CCA')
with SM.sl(['chi_ang_trans'],ds,1,mode):
	chi_ang_trans=np.cos(chi_ang[0])
	chi_ang_trans=np.vstack((chi_ang_trans,np.sin(chi_ang[0])))
	for i in chi_ang[1:]:
		chi_ang_trans=np.vstack((chi_ang_trans,np.cos(i)))
		chi_ang_trans=np.vstack((chi_ang_trans,np.sin(i)))
with SM.sl(['CCA_rslt'],ds,1,mode):
	CCA_rslt=CCA.cc(D.T,chi_ang_trans.T)
CCA_eigva=CCA_rslt.rx2('cor')
CCA_eigveX=CCA_rslt.rx2('xcoef') #Distance
CCA_eigveY2=np.array(CCA_rslt.rx2('ycoef')).T #Chi_ang in terms of sin and cos
CCA_eigveY=[]
for i in CCA_eigveY2:
	CCA_eigveY.append(sqrt(i[::2]**2+i[1::2]**2))
CCA_eigveY=np.vstack(CCA_eigveY)
CCA_eigveY=CCA_eigveY.T

#-----------Secondary Calculations----------------
#----------------Plotting-------------------------
print 'Plotting my results..    ('+tim()+')'
#-------------------------------
#	Scatter and Projections
#-------------------------------
##{{{
for i in xrange(4):
	# http://matplotlib.org/examples/pylab_examples/scatter_hist.html
#=====================
#   Setting plots up
#=====================
	##{{{
	#----------Setting up Dimensions
	width1,width2,width3=0.12,0.58,0.06
	hgap1,hgap2=0.005,0.11
	left1= 0.05
	left2=left1+width1+hgap1
	left3=left2+width2+hgap2
	height1,height2=0.62,0.15
	vgap1=0.008
	height3=height1+vgap1+height2
	bottom1=0.110
	bottom2=bottom1+height1+vgap1
	bottom3=bottom1
	rect_v1v2= [left2, bottom1, width2, height1]
	rect_v1 = [left1, bottom1, width1, height1]
	rect_v2 = [left2, bottom2, width2, height2]
	rect_cb=[left3,bottom3,width3,height3]
	fig=vars()['fig'+str(i)]=plt.figure(i,figsize=(11.69, 8.27),dpi=DPI)
	##}}}
#=====================
#	Plotting Vectors
#=====================
	##{{{
	ax_v1v2 = plt.axes(rect_v1v2)
	ax_v1 = plt.axes(rect_v1,sharey=ax_v1v2)
	ax_v2 = plt.axes(rect_v2,sharex=ax_v1v2)
	# Labels and titles
	ax_v1.get_yaxis().set_visible(False) #make ticks invisible
	ax_v2.get_xaxis().set_visible(False) 
	plt.setp(ax_v1.xaxis.get_ticklabels(),size=7,rotation=90) #Rotate ticks and change font size
	plt.setp(ax_v2.yaxis.get_ticklabels(),size=7)
	ax_v1.set_title('Displacement from ave',position=(-0.1,0.5),transform=ax_v1.transAxes,rotation='vertical',fontsize='9')
	ax_v2.set_title('Chi_Angle',fontsize='9')
	ax_v2.set_ylabel('Contribution to vector',fontsize='7')
	nn=np.array
	#-----------------------------------------------------------------------------
	# Side graphs
	v1=abs(nn(CCA_eigveX)[:,i]/max(nn(CCA_eigveX)[:,i]))	#Distances normalized into percentages
	v1=v1/sum(v1)
	bar11=ax_v1.barh(arange(len(v1)),v1,align='center',edgecolor='#ff5252',color='#ffcdcd')
	v2=abs(nn(CCA_eigveY)[:,i]/max(nn(CCA_eigveY)[:,i]))	#chi_ang normalized
	v2=v2/sum(v2)
	bar21=ax_v2.bar(arange(len(v2)),v2,align='center',edgecolor='#ff5252',color='#ffcdcd')
	ax_v1.set_xlim(ax_v1.get_xlim()[::-1]) #Flip bar graph
	#-----------------------------------------------------------------------------
	# Main Graphs
	X,Y=np.meshgrid(v2,v1)
	Z=(X*Y)
	surf=ax_v1v2.pcolor(Z,cmap='Reds',rasterized=True)
	ax_v1v2.autoscale(enable=True, axis='both', tight=True)
	Ylabel=struc.ca_names #distance
	Xlabel=chi_resid #chi_angle
	# Set xaxis ticks
	ax_v1v2.set_xticks(arange(len(Xlabel)))
	ax_v1v2.set_xticklabels(Xlabel,size=5)
	plt.setp(ax_v1v2.xaxis.get_ticklabels(),va='top',ha='left',rotation=90) #Rotates
	# Set yaxis ticks
	ax_v1v2.yaxis.tick_right()
	ax_v1v2.set_yticks(arange(len(Ylabel)))
	ax_v1v2.set_yticklabels(Ylabel,size=4)
	plt.setp(ax_v1v2.yaxis.get_ticklabels(),va='bottom',ha='left')
	#ax_v1v2.grid(True)
	#-----------------------------------------------------------------------------
	# Colorbar and creating color bar axes
	ax_cm=plt.axes(rect_cb)
	cb=plt.colorbar(surf,cax=ax_cm)
	plt.setp(ax_cm.get_xticklabels(),size=6)
	ax_cm.yaxis.set_ticks(np.linspace(ax_cm.yaxis.get_ticklocs()[0],ax_cm.yaxis.get_ticklocs()[-1],3))
	ax_cm.yaxis.set_ticklabels(['low','med','high'])
	fig.text(0.1,0.96,os.getcwd(),fontsize=9)
	fig.text(0.1,0.92,'correlation coefficient = {0:0.2f}'.format(CCA_eigva[i]),fontsize=9)
	MPL.show()
	##}}}
##}}}
#---------------------
#	Plotting Scatter
#---------------------
##{{{
fig8=plt.figure(8,figsize=(11.69,8.27),dpi=DPI)
for i in xrange(4):
	ax_cluster=fig8.add_subplot(2,2,i+1)
	print 'Projecting data onto eigenvectors {0:0.0f}...    ('.format(i+1)+tim()+')'
	proj_dis=np.dot(D.T,nn(CCA_eigveX)[:,i]) #distance
	proj_chi=np.dot(chi_ang.T,nn(CCA_eigveY)[:,i]) #chi_ang
	# Plotting clusters
	lin=ax_cluster.scatter(proj_chi,proj_dis,c=arange(len(proj_dis)),cmap='gist_rainbow',marker='o',lw=0,s=2, rasterized=True)
	ax_cluster.set_ylabel('Distance proj'+str(i+1),fontsize=10)
	ax_cluster.set_xlabel('Chi Angle proj'+str(i+1),fontsize=10)
	cb=plt.colorbar(lin,ax=ax_cluster)
	cb.set_label('Frames',fontsize=8)
	cb.set_ticks([0,1])
	cb.set_ticklabels(['0','.'])
fig8.text(0.1,0.96,os.getcwd(),fontsize=9)
MPL.show()
##}}}

#----------------------
#	Correlation bar
#----------------------
##{{{
fig7=plt.figure(7,figsize=(11.69, 8.27),dpi=DPI)
ax17=fig7.add_subplot(111)
b117=ax17.bar(arange(len(CCA_eigva)),CCA_eigva,width=0.8,color='#9f9f9f')
ax17.set_xlabel('Eigenvector number')
ax17.set_ylabel('Correlation coefficient value')
ax17.set_title('Canonical Correlation Analysis')
MPL.show()
##}}}

##}}}



