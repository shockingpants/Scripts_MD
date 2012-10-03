#!/opt/local/bin/python2.7
# Calc_Delta.py

from numpy import *
import numpy as np
import scipy
import scipy.stats as ss
import os
import sys
import pprint
import matplotlib.pyplot as plt
import matplotlib
import MDAnalysis
from strtup import *
import math
from mpl_toolkits.axes_grid1 import host_subplot
from mpl_toolkits.mplot3d import Axes3D 
from matplotlib.backends.backend_pdf import PdfPages

##################################
####      Extracting data
##################################
##{{{
#===========Load data files and setting class
# The names are Complex Receptor Ligand
try:
	sol=sys.argv[1]
except Exception as err:
	print "(JON)Warning:", err, " Using Default instead."
	sol="GB" #Set which solvent approximation we are using
class energy:
	def __init__(self):
		self.TYPES=["Bond","Angle","Dihedral","OnetoFourvdw","VDWAALS","EEL","ESURF","E"+sol]
		self.HEADERS=["Complex","Receptor","Ligand"]
		for i in range(len(self.TYPES)):
			vars(self)[self.TYPES[i]]=genfromtxt(self.TYPES[i]+'.dat', names=True, dtype=('float', 'float', 'float'))
			if not 'frames' in vars(self):
				self.frames=arange(len(vars(self)[self.TYPES[i]]))
		self.be()
	def be(self):
		# Binding energy
		for i in range(len(self.TYPES)):
			vars(self)[self.TYPES[i]+"_BE"]=vars(self)[self.TYPES[i]][self.HEADERS[0]]-vars(self)[self.TYPES[i]][self.HEADERS[1]]-\
			vars(self)[self.TYPES[i]][self.HEADERS[2]]
			if i == 0:
				self.BE=vars(self)[self.TYPES[i]+"_BE"]
			else: 
				self.BE=self.BE+vars(self)[self.TYPES[i]+"_BE"]
#===========Saving to file
A=energy() #Instantiate all relavant energy
headers=['Index','BE']
headers.extend(A.TYPES)
fid=open('BE_Timeseries','w')
fid.write(' '.join(headers))
fid.write("\n")
for i in range(len(A.frames)):
	# print out line by line
	fid.write('{0:d} {1:5f} {2:5f} {3:5f} {4:5f} {5:5f} {6:5f} {7:5f} {8:5f} {9:5f}\n'.format(i,A.BE[i],*[vars(A)[A.TYPES[j]+"_BE"][i] for j in range(8)]))
fid.close()
##}}}
##################################
####	   Plotting data
##################################
# To allow both shell run and  interactive execfile
try:
	closeall()
except NameError:
	pass
#======== Set up pages
pp = PdfPages('Binding_Energy.pdf')
DPI=200
#plt.ion()
# We 4 main energies and we will plot that along with their correlation
EOI=["VDWAALS_BE", "EEL_BE", "ESURF_BE", "E"+sol+"_BE"]
n=30 #skip
#======== Page 1
#{{{
#=================
#  Binding Energy
#=================
fig3=plt.figure(3,figsize=(11.69, 8.27),dpi=DPI)
ax1=fig3.add_subplot(3,1,1)
ax1.set_title('Binding Energy')
lin1, =ax1.plot(A.frames,vars(A)['BE'])
lin2, =ax1.plot(A.frames,mean(A.BE)*np.ones(shape(A.frames)),'-r')
lin3, =ax1.plot(A.frames,(mean(A.BE)+std(A.BE))*np.ones(shape(A.frames)),'--r')
lin4, =ax1.plot(A.frames,(mean(A.BE)-std(A.BE))*np.ones(shape(A.frames)),'--r')
ax1.legend([lin2,lin3],['Mean', '1 STD'])
ax1.get_legend().get_frame().set_alpha(0.3)
ax1.set_ylabel('kcal mol$^{-1}$')
ax1.get_xaxis().set_visible(False)
ax1.text(0,1.2,os.getcwd(),fontsize=9,transform=ax1.transAxes) # FOLDER NAME
#===========================
#  MVA and other smoothing
#===========================
#===== MVA
window=30
ax2=fig3.add_subplot(312,sharey=ax1)
lin1, =ax2.plot(A.frames,vars(A)['BE'],color='yellow')
lin2, =ax2.plot(A.frames[::n][mva.sma(vars(A)['BE'][::n],window)[0]],mva.sma(vars(A)['BE'][::n],window)[1])
lin3, =ax2.plot(A.frames[::n][mva.wma(vars(A)['BE'][::n],window)[0]],mva.wma(vars(A)['BE'][::n],window)[1])
#===== FFT
lin4, =ax2.plot(*fft_mfreq(A.BE,A.frames,high=0.005,low=0))
MVA=mva.sma(vars(A)['BE'][::n],window)[1]
lin5, =ax2.plot(A.frames[::n][mva.sma(vars(A)['BE'][::n],window)[0]],mva.sma(vars(A)['BE'][::n],window)[1]+std(A.BE),'--g')
lin6, =ax2.plot(A.frames[::n][mva.sma(vars(A)['BE'][::n],window)[0]],mva.sma(vars(A)['BE'][::n],window)[1]-std(A.BE),'--g')
#===== PLots
leg2=ax2.legend([lin1,lin2,lin3,lin4],['BE','SMA','WMA','FFT'])
ax2.set_title('Moving Averages')
ax2.set_ylabel('kcal mol$^{-1}$')
#ax2.set_xlabel('Frame')
#ax2.get_xaxis().set_visible(False)
plt.setp(ax2.get_xticklabels(),fontsize=10)
ax2.get_legend().get_frame().set_alpha(0.3)
ax2.text(0.03,0.95,r'$\mu$={0:4.4f}, $\sigma$={1:4.4f}, skip={2:1d}, window={3:1d}, Total_time={4:0.0f}ns'\
.format(mean(A.BE), std(A.BE), n, window, len(A.frames)*0.002), fontsize=10,transform=ax2.transAxes,verticalalignment='top')
#===========================
#	Auto Correlation
#===========================
ax3=fig3.add_subplot(313)
ax3.plot(*autocorr(A.BE))
ax3.set_title('Autocorrelation')
plt.setp(ax3.get_xticklabels(),fontsize=10)
pp.savefig(format='pdf')
##}}}
#======== Page 2
##{{{
#== Correlations
for i in range(4):
	vars()['coor'+str(i)]=np.corrcoef(A.BE,vars(A)[EOI[i]])[0,1]
	vars()['MEAN'+str(i)]=mean(vars(A)[EOI[i]])
	vars()['STD'+str(i)]=std(vars(A)[EOI[i]])

fig1=plt.figure(1,figsize=(11.69, 8.27),dpi=DPI)
#plt.ion()
#plt.show()
#== Create Subplots
ax1=fig1.add_subplot(4,1,1)
ax1a=ax1.twinx()
ax2=fig1.add_subplot(4,1,2) 
ax2a=ax2.twinx()
ax3=fig1.add_subplot(4,1,3)
ax3a=ax3.twinx()
ax4=fig1.add_subplot(4,1,4)
ax4a=ax4.twinx()
#===Plot
lin1a, = ax1a.plot(A.frames[::n],vars(A)['BE'][::n],'r')
plt.setp(lin1a,alpha=0.7)
lin1, = ax1.plot(A.frames[::n],vars(A)[EOI[0]][::n],'b')
lin2a, = ax2a.plot(A.frames[::n],vars(A)['BE'][::n],'r')
plt.setp(lin2a,alpha=0.7)
lin2, = ax2.plot(A.frames[::n],vars(A)[EOI[1]][::n],'b')
lin3a, = ax3a.plot(A.frames[::n],vars(A)['BE'][::n],'r')
plt.setp(lin3a,alpha=0.7)
lin3, = ax3.plot(A.frames[::n],vars(A)[EOI[2]][::n],'b')
lin4a, = ax4a.plot(A.frames[::n],vars(A)['BE'][::n],'r')
plt.setp(lin4a,alpha=0.7)
lin4, = ax4.plot(A.frames[::n],-vars(A)[EOI[3]][::n],'b')
#== Insert Labels
ax1.set_ylabel("kcal mol-1")
ax2.set_ylabel("kcal mol-1")
ax3.set_ylabel("kcal mol-1")
ax4.set_ylabel("kcal mol-1")
#== Insert Corr
for i in range(4):
	COOR=vars()['coor'+str(i)]
	MEAN=vars()['MEAN'+str(i)]
	STD=vars()['STD'+str(i)]
	vars()['ax'+str(i+1)].text(0.95,0.05,r'coor={0:0.2f},$\mu$={1:0.2f},$\sigma$={2:0.2f}'.format(COOR,MEAN,STD),\
	horizontalalignment='right',verticalalignment='bottom',transform=vars()['ax'+str(i+1)].transAxes)
#== Remove ticks
ax1.get_xaxis().set_visible(False)
ax2.get_xaxis().set_visible(False)
ax3.get_xaxis().set_visible(False)
#===Insert Legend
#leg1=ax1a.legend([lin1, lin1a],[EOI[0], "BE"])
#leg2=ax2a.legend([lin2, lin2a],[EOI[1], "BE"])
#leg3=ax3a.legend([lin3, lin3a],[EOI[2], "BE"])
#leg4=ax4a.legend([lin4, lin4a],["(-ve)"+EOI[3], "BE"])
#leg1.get_frame().set_alpha(0.3)
#leg2.get_frame().set_alpha(0.3)
#leg3.get_frame().set_alpha(0.3)
#leg4.get_frame().set_alpha(0.3)
#===Titles and stuff
ax1.set_title(EOI[0],fontsize=12)
ax2.set_title(EOI[1],fontsize=12)
ax3.set_title(EOI[2],fontsize=12)
ax4.set_title("(-ve)"+EOI[3],fontsize=12)
yran=ax1.get_ylim()
ax1.text(0,1.2,os.getcwd(),fontsize=9,transform=ax1.transAxes) # FOLDER NAME
pp.savefig(format='pdf')
##}}}
#======== Page 3 (Histogram of energy)
##{{{
fig2=plt.figure(2,figsize=(11.69, 8.27),dpi=DPI)
ax1=fig2.add_subplot(111)
# the histogram of the data
n, bins, patches = ax1.hist(A.BE, 150, normed=1, facecolor='green', alpha=0.75)#======= Page 4
ax1.set_xlabel(r'BE kcal mol$^{-1}$')
ax1.set_ylabel('Frequency')
def func(x,u,s):
	return (1./s/sqrt(2*pi))*exp(-(x-u)**2./(2.*s**2.))
bins=diff(bins)/2+bins[0:-1]
x, fit_y, GOF, fit=curve_fit(func,bins,n,-60,15)
lin2, = ax1.plot(x,fit_y,'r--')
ax1.text(0.8,0.95,'$\mu$ = {0:0.3f}\n$\sigma$ = {1:0.3f}'.format(fit.param[0],fit.param[1]),horizontalalignment='left',verticalalignment='top',transform=ax1.transAxes)
ax1.text(0,1.05,os.getcwd(),fontsize=9,transform=ax1.transAxes) # FOLDERNAME
pp.savefig(format='pdf')
##}}}
pp.close()





 

