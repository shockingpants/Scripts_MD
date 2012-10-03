!/opt/local/bin/python2.7
# struc_ana.py
# Calculates structural Related Quantities in a protein
# Note, Pretty much all the analysis depends on R and bio3d.
# Can implement in MDAnalysis... but...
#=============== MODULES ======================
##{{{
from numpy import *
import numpy as np
import scipy
import scipy.stats as ss
import os
import sys
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


##}}}
#==============================================

########################################
####        Script Startup
########################################
b3d=importr('bio3d')
mode=[] # 0 --> run only 00 --> run or load  1 --> run and save  2 --> load
filename='dat0.db'
mode,filename = SM.saveload(mode, filename) # Run this for checks.
if mode == 2:
	ds=SM.get_set(filename)
start_time=time.time()
#mmode decides if the small sections will be loaded or not

########################################
####	Parameters
########################################
##{{{ Param
# Contains Mainly File Locations
print 'Loading parameters...'
if mode == 1 or mode == 0:
	param={}
	param['dssp']="/opt/local/bin/"
	param['trj_folder']='../../' # Where the trajectories are located relative to script
elif mode == 2:
	param=ds.param

for i in iter(param):
	vars()[i]=param[i]

##}}}
########################################
####    Loading files and extract
########################################
##{{{ Trajectories
print 'Loading files and trajectories...'
#==== Read a trajectory file
trj=b3d.read_ncdf(trj_folder+"prod_CA.nc")
#=== Read a PDB file with the CA from the minimized structure
pdb=b3d.read_pdb(trj_folder+"prod_CA.pdb.1")
ca_inds=b3d.atom_select(pdb, "calpha")
pdb2=b3d.read_pdb(trj_folder+"prot.pdb") # Crystal Structure
ca_inds2=b3d.atom_select(pdb2, "calpha").rx2('atom') # Crystal Structure with calpha only
ca_seq=b3d.aa321(pdb2.rx2('atom').rx(ca_inds2,'resid'))
ca_resno=pdb2.rx2('atom').rx(ca_inds2,'resno')
ca_names=R.paste(ca_seq,ca_resno) # x-axis labels
ca_chain=pdb2.rx2('atom').rx(ca_inds2,'chain')
ca_len=R.length(ca_inds2)
ca_breaks=R.list()
##}}}
########################################
####    Analysis
########################################
##{{{ Analysis
print 'Performing Analysis...  ('+tim()+')'
#=======================================
#   Calculating Secondary Structure
#=======================================
##{{{SSE
print 'Assigning SSE...  ('+tim()+')'
mmode=SM.msaveload(0,mode)
#----------------------------------------------
if mmode == 0:
	sse=b3d.dssp(pdb2, exepath = "/opt/local/bin/", resno=False)
	# Assign helixes
	if not R('is.null')(sse.rx2('helix').rx2('start'))[0]: #returns True
		sse.rx2('helix').rx2('start').ro=sse.rx2('helix').rx2('start').ro-0.5
		sse.rx2('helix').rx2('end').ro=sse.rx2('helix').rx2('end').ro+0.5
	# Assign sheets
	if not R('is.null')(sse.rx2('sheet').rx2('start'))[0]: #returns True
		sse.rx2('sheet').rx2('start').ro=sse.rx2('sheet').rx2('start').ro-0.5
		sse.rx2('sheet').rx2('end').ro=sse.rx2('sheet').rx2('end').ro+0.5

elif mmode==2:
	print 'Loading...'
	sse=ds.sse

##}}}SSE
#=======================================
#          RMSD, RMSF and Fitting 
#=======================================
##{{{ RMSD,RMSD and FITTING 
print 'Calculating RMSD, RMSF and fitting...  ('+tim()+')'
mmode=SM.msaveload(0,mode)
#------------------------------------------------------------
if mmode==0:
	#====== Fit trajectory onto pdb with CA
	print 'Fitting...   ('+tim()+')'
	# returns a <matrix> in R
	xyz=b3d.fit_xyz(pdb.rx2('xyz'), trj, fixed_inds=ca_inds.rx2('xyz'), mobile_inds=ca_inds.rx2('xyz')) #This will be a VERY long step. May be faster in... MDAnalysis?
	#====== Getting RMSD
	print 'RMSDs...  ('+tim()+')'
	# all return <numeric> in R
	ave=R.apply(xyz,2,R.mean) # Gets the average Structure
	r_ave=b3d.rmsd(a=ave, b=xyz) # Slow # Gets RMSD from ave structure over time
	r_first=b3d.rmsd(a=xyz.rx(1,True), b=xyz) # Slow # RMSD from first structure
	r_crystal=b3d.rmsd(a=R.c(R('as.numeric')(R.t(pdb2.rx2('atom').rx(ca_inds2,R.c('x','y','z'))))), b=xyz, fit=True) # Gets RMSD from crystal Structure
	#====== Getting RMSF
	print 'RMSFs...   ('+tim()+')'
	fluc=b3d.rmsf(xyz)
	

elif mmode==2:
	print 'Loading...'
	xyz=ds.xyz
	ave=ds.ave
	r_ave=ds.r_ave
	r_first=ds.r_first
	r_crystal=ds.r_crystal
	fluc=ds.fluc
##}}}
#=======================================
#         DCCM  and  XYZ Corr 
#=======================================
##{{{
print 'Calculating Correlation Matrix based on each rmsd (DCCM)...  ('+tim()+')'
mmode=SM.msaveload(0,mode)
#-------------------------------------------------
if mmode==0:
	cij=b3d.dccm(xyz) # numpy's implementation is faster.... 
elif mmode==2:
	print 'Loading...'
	cij=ds.cij 

print 'Calculating Correlation Matrix based on each xyz...  ('+tim()+')'
mmode=SM.msaveload(0,mode)
#-------------------------------------------------
if mmode==0:
	XYZ=np.array(xyz)# Rows=degrees of freedom(3xno. of CA) # Columns=time 
	AVE=np.array(ave)
	A,B=np.meshgrid(AVE,xrange(shape(XYZ)[0]))
	net_XYZ=XYZ-A # This is the individual x,y and z distances from the mean structure
	corr_xyz=np.corrcoef(net_XYZ.T)
elif mmode==2:
	print 'Loading...'
	XYZ=ds.XYZ
	AVE=ds.AVE
	corr_xyz=ds.corr_xyz
##}}}
#=======================================
#               PCA
#=======================================
##{{{
print 'Calculating PCA...  ('+tim()+')'
mmode=SM.msaveload(0,mode)
#-------------------------------------------------
if mmode==0:
	pca_raw=b3d.pca_xyz(xyz)
	eigenvectors=pca_raw.rx2('U')
	projected_data=pca_raw.rx2('z')
	scores=pca_raw.rx2('L') # This is also the eigenvalues!
elif mmode==2:
	print 'Loading...'
	pca_raw=ds.pca_raw
	eigenvectors=ds.eigenvectors
	projected_data=ds.projected_data
	scores=ds.scores

##}}}
#=======================================
#            Saving Data
#=======================================
##{{{Save
if mode == 1:
	print 'Saving data...   ('+tim()+')'
	t=time.time()
	ds=SM.add_set()
	ds.save_all('sse', 'xyz', 'ave', 'r_ave', 'r_first', 'r_crystal', 'cij', 'fluc','XYZ', 'AVE', 'corr_xyz', 'pca_raw', 'eigenvectors', 'projected_data', 'scores', sse, xyz, ave, r_ave, r_first, r_crystal, cij, fluc, XYZ, AVE, corr_xyz, pca_raw, eigenvectors, projected_data, scores, param=param, mode=1)
	print 'Finished saving in {0:0.1f}s.  ('.format(time.time()-t)+tim()+')'
else:
	print 'All variables have been loaded. Nothing is saved. ('+tim()+')'
##}}}

#Time Check
print '{0:0.1f}s has elapsed.'.format(time.time()-start_time)
##}}}
########################################
####    Plotting
########################################
##{{{Plotting
print 'Plotting Graphs...    ('+tim()+')'
MPL=MM.mpl([],'struc_ana.pdf')
DPI=150
#=====================
#   Load Colomaps
#=====================
##{{{
cpool = [ '#0308d3','#0308d3','#4b47fe','#4b47fe', '#9393f4','#9393f4','#dbdbfb', '#dbdbfb','#ffffff','#ffffff','#ffffff','#ffffff','#ffe1e1','#ffe1e1','#fd8b8b','#fd8b8b','#f95454', '#f95454','#ff0000','#ff0000']
my_cmap = matplotlib.colors.ListedColormap(cpool[0:], 'indexed')
##}}}
#=====================
#    Figure 1 
#=====================
##{{{ DCCM
print 'Plotting correlation between distances from avg struct...    ('+tim()+')'
fig1=plt.figure(1,figsize=(11.69, 8.27),dpi=DPI)
ax11=fig1.add_subplot(111)
#ele111=ax11.pcolor(np.array(cij),cmap=my_cmap,vmax=1,vmin=-1,rasterized=True)
ele111=ax11.contourf(np.array(cij),arange(-1,1.2,0.2),cmap=my_cmap,rasterized=True)
ax11.set_title('Correlation of distances from avg structure over time')
ax11.set_xticks(arange(len(np.array(ca_names))))
ax11.set_xticklabels(np.array(ca_names),fontsize=5,ha='center')
ax11.set_yticks(arange(len(np.array(ca_names))))
ax11.set_yticklabels(np.array(ca_names),fontsize=5,ha='right')
fig1.autofmt_xdate(ha='left',rotation=90)
ax11.set_xlabel('Residues')
ax11.grid(True, which='major',ls='-',color='w')
plt.colorbar(ele111)
fig1.text(0.1,0.96,os.getcwd(),fontsize=9)
MPL.show()
##}}}
#=====================
#    Figure 2
#=====================
##{{{XYZ Corr
print 'Plotting Correlation between XYZ...    ('+tim()+')'
fig2=plt.figure(2,figsize=(11.69, 8.27),dpi=DPI)
ax21=fig2.add_subplot(111)
ele211=ax21.pcolor(corr_xyz,cmap=my_cmap,vmax=1,vmin=-1,rasterized=True)
ax21.set_xticks(arange(shape(corr_xyz)[0])[::3])
ax21.set_xticklabels(np.array(ca_names),fontsize=5,ha='right')
ax21.set_yticks(arange(shape(corr_xyz)[0])[::3])
ax21.set_yticklabels(np.array(ca_names),fontsize=5,ha='right')
fig2.autofmt_xdate(ha='left',rotation=90)
ax21.set_xlabel('Residues')
ax21.grid(True, which='major',ls='-',color='w')
plt.colorbar(ele211)
fig2.text(0.1,0.96,os.getcwd(),fontsize=9)
MPL.show()
##}}}
#=====================
#    Figure 3
#=====================
##{{{RMSD
print 'Plotting Various RMSD over time...   ('+tim()+')'
fig3=plt.figure(3,figsize=(11.69, 8.27),dpi=DPI)
RMSD_dic={'r_ave':r_ave,'r_crystal':r_crystal,'r_first':r_first}
for ind,(i,RMSD) in enumerate(RMSD_dic.iteritems()):
	RMSD=np.array(RMSD)
	frames=arange(len(RMSD))
	ax=fig3.add_subplot(len(RMSD_dic.keys()),1,ind+1)
	ax.set_ylabel('RMSD (A)')
	ax.set_title(i)
	lin,=ax.plot(frames[::10],RMSD[::10])
	lin2,=ax.plot(frames[::10][mva.sma(RMSD[::10],100)[0]],mva.sma(RMSD[::10],100)[1]) #Simple Average over 10%
	lin3,=ax.plot(frames[::10][mva.wma(RMSD[::10],100)[0]],mva.wma(RMSD[::10],100)[1]) #Weighted Average over 10%
	lin4,=ax.plot(*fft_mfreq(RMSD,frames,high=0.001,low=0)) #FFT Fit
	lin5,=ax.plot(frames[::10][mva.sma(RMSD[::10],100)[0]],mva.sma(RMSD[::10],100)[1]+std(RMSD),'--g') # SMA+std
	lin6,=ax.plot(frames[::10][mva.sma(RMSD[::10],100)[0]],mva.sma(RMSD[::10],100)[1]-std(RMSD),'--g') # SMA-std
	ax.set_xlim(0,max(frames))
	leg=ax.legend((lin,lin2,lin3,lin4),(i,'SMA','WMA','FFT'),prop={'size':9})
	leg.get_frame().set_alpha(0.3) # Setting Transparency
	x=ax.get_xaxis()
	x.set_visible(False)
	ax.text(0.03,0.95,r'$\mu$={0:4.4f}, $\sigma$={1:4.4f}, skip={2:1d}, window={3:1d}, Total_time={4:0.0f}ns'\
	.format(mean(RMSD), std(RMSD), 10, 100, max(frames)*0.002), fontsize=10,transform=ax.transAxes,verticalalignment='top')
	vars()['ax'+'3'+str(ind+1)]=ax
	vars()['leg'+'3'+str(ind+1)]=leg
# Sets property for last axes
x.set_visible(True)
ax.set_xlabel('Frames')
fig3.text(0.1,0.96,os.getcwd(),fontsize=9)
MPL.show()
##}}}
#=====================
#    Figure 4
#=====================
##{{{ RMSF
print 'Plotting RMSF over residues...  ('+tim()+')'
fig4=plt.figure(4,figsize=(11.69, 8.27),dpi=DPI)
width = 1
left = range(len(np.array(fluc)))
ax41=fig4.add_subplot(111)
#bar=ax41.bar(left, np.array(fluc), width, color="#ffe0bd",edgecolor="#ffc27a")
lin1,=ax41.plot(np.array(fluc),color='#000000')
ax41.set_xticks(arange(len(np.array(fluc)))+0.5)
ax41.set_xticklabels(np.array(ca_names),fontsize=7,ha='right')
fig4.autofmt_xdate(ha='left',rotation=90)
ax41.set_xlabel('Residues')
ax41.set_ylabel('RMSF (A)')
ax41.set_title('Root Mean Squared Fluctuations')
# Applying Secondary Structures
helix=[(sse.rx2('helix').rx2('start')[i],sse.rx2('helix').rx2('end')[i]) for i in xrange(len(sse.rx2('helix').rx2('start')))]
sheet=[(sse.rx2('sheet').rx2('start')[i],sse.rx2('sheet').rx2('end')[i]) for i in xrange(len(sse.rx2('sheet').rx2('start')))]
for i in helix: 
	if i == helix[0]:
		plt.axvspan(i[0],i[1],color='r',alpha=0.5,label='helix')
	else:
		plt.axvspan(i[0],i[1],color='r',alpha=0.5)

for i in sheet:
	if i == sheet[0]:
		plt.axvspan(i[0],i[1],color='#fbfe00',alpha=0.5,label='sheet')
	else:
		plt.axvspan(i[0],i[1],color='#fbfe00',alpha=0.5)

plt.legend()
fig3.text(0.1,0.96,os.getcwd(),fontsize=9)
MPL.show()
##}}}
#=====================
#    Figure 5,6,7,8
#=====================
##{{{ PCA
#===== Figure 5
print 'Plotting data projected onto various PCA...  ('+tim()+')'
fig5=plt.figure(5,figsize=(11.69, 8.27),dpi=DPI)
PCA_proj=np.array(projected_data)[:,0:6]
PCA_proj=PCA_proj.T
for i in xrange(6):
	ax=fig5.add_subplot(3,2,i+1)
	ax.plot(PCA_proj[i])
	ax.set_xlim(0,len(PCA_proj[i]))
	ax.set_title('Projected onto PCA'+str(i+1),fontsize=9)
	ax.set_xlabel('Frames')
	ax.set_ylabel('Fluctuations')
	if not (i==4 or i==5):
		ax.get_xaxis().set_visible(False)

fig5.text(0.1,0.96,os.getcwd(),fontsize=8)
fig5.text(0.5,0.92,'Data Projection onto PCAs',fontsize=12,ha='center')
MPL.show()

#===== Figure 6
print 'Plotting contribution by components in PCA...  ('+tim()+')'
fig6=plt.figure(6,figsize=(11.69, 8.27),dpi=DPI)
E_vectors=np.array(eigenvectors)
E_vectors=E_vectors.T
E_sum=sum(abs(E_vectors),axis=1)
temp,E_sum2=meshgrid(E_sum,E_sum)
EV_norm=abs(E_vectors)/E_sum2 # Normalizing
for i in xrange(3):
	ax=fig6.add_subplot(3,1,i+1)
	lin, =ax.plot(EV_norm[i,:],alpha=0.1)
	lin2, =ax.plot(EV_norm[i,:],'r.')
	ax.set_xlim(0,shape(EV_norm)[0])
	ax.set_ylim(0)
	ax.set_title('Coordinate contribution of PCA'+str(i+1),fontsize=9)
	ax.set_ylabel('Contribution Proportion')
	ax.set_xlabel('Residues XYZ')
	ax.set_xticks(arange(shape(EV_norm)[0])[::3])
	ax.set_xticklabels(np.array(ca_names),fontsize=7,ha='right')
	fig6.autofmt_xdate(ha='left',rotation=90)
	for i in (range(shape(EV_norm)[0])[::6]):
		plt.axvspan(i-0.4,i+2.4,color='#fea100',alpha=0.2)

fig6.text(0.1,0.96,os.getcwd(),fontsize=9)
MPL.show()
#==== Figure 7
fig7=plt.figure(7,figsize=(11.69, 8.27),dpi=DPI)
E_vectors=np.array(eigenvectors)
E_vectors=E_vectors.T
E_sum=sum(abs(E_vectors),axis=1)
temp,E_sum2=meshgrid(E_sum,E_sum)
EV_norm=abs(E_vectors)/E_sum2 # Normalizing
for i in xrange(3,6):
	ax=fig7.add_subplot(3,1,i-3+1)
	lin, =ax.plot(EV_norm[i,:],alpha=0.1)
	lin2, =ax.plot(EV_norm[i,:],'r.')
	ax.set_xlim(0,shape(EV_norm)[0])
	ax.set_ylim(0)
	ax.set_title('Coordinate contribution of PCA'+str(i+1),fontsize=9)
	ax.set_ylabel('Contribution Proportion')
	ax.set_xlabel('Residues XYZ')
	ax.set_xticks(arange(shape(EV_norm)[0])[::3])
	ax.set_xticklabels(np.array(ca_names),fontsize=7,ha='right')
	fig7.autofmt_xdate(ha='left',rotation=90)
	for i in (range(shape(EV_norm)[0])[::6]):
		plt.axvspan(i-0.4,i+2.4,color='#fea100',alpha=0.2)

fig6.text(0.1,0.96,os.getcwd(),fontsize=9)
MPL.show()

#==== Figure 8
print 'Plotting scatter plot between PCA 1--2, 1--3, 2--3...  ('+tim()+')'
fig8=plt.figure(8,figsize=(11.69, 8.27),dpi=DPI)
# Plotting clusters
order={0:(0,1),1:(0,2),2:(1,2)}
frame=arange(len(PCA_proj[1]))
for i,(a,b) in order.iteritems(): 
	ax=fig8.add_subplot(2,2,i+1)
	lin=ax.scatter(PCA_proj[a],PCA_proj[b],c=frame,cmap='gist_rainbow',marker='o',lw=0,s=2, rasterized=True)
	ax.set_xlabel('PCA'+str(a+1),fontsize=10)
	ax.set_ylabel('PCA'+str(b+1),fontsize=10)
	cb=plt.colorbar(lin)
	cb.set_label('Frames',fontsize=8)
	cb.set_ticks([0,1])
	cb.set_ticklabels(['0','.'])
# Plotting PCA amount
scoresum=sum(scores)
scoreprop=scores/scoresum # Proportion
ax84=fig8.add_subplot(2,2,4)
x_ax=arange(10)+1
lin841,=ax84.plot(x_ax,scoreprop[0:10],alpha=0.2)
lin842,=ax84.plot(x_ax,scoreprop[0:10],'dr')
ax84.set_xlabel('PCA components',fontsize=10)
ax84.set_ylabel('Proportion',fontsize=10)
for ind,i in enumerate(scoreprop[0:8]):
	scorecum=sum(scoreprop[0:ind+1])*100
	ax84.annotate('{0:0.0f}%'.format(scorecum),xy=(x_ax[ind],scoreprop[ind]), xytext=(x_ax[ind]+0.2, scoreprop[ind]),fontsize=9)

fig8.text(0.1,0.96,os.getcwd(),fontsize=9)
MPL.show()
##}}}
##}}}
print '{0:0.1f}s has elapsed'.format(time.time()-start_time)

##}}}


