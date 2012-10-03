#!/opt/local/bin/python2.7
# plot_ala.py
# import numpy as np
from numpy import *
# Make all matlib functions accessible at the top level via M.func()
import numpy.matlib as M
# Make some matlib functions accessible directly at the top level via, e.g. rand(3,3)
from numpy.matlib import rand,zeros,ones,empty,eye
import os
import sys
import pprint
import matplotlib.pyplot as plt

################
### Import Files
################
# Headers: resn Tot_energy STD STDE_mean
mut_ddG=genfromtxt('mut_ddG.dat',skip_header=0, names=True, delimiter=" ", dtype=('|S5, float'))
# Title
if len(sys.argv)==2:
	title_1 = sys.argv[1]
	title_1='('+title_1+')'
else:
	title_1=''

###########################
###      Plot Graph
###########################
##{{{
#==========
#   Plot
#==========
ener=mut_ddG['ddG']
#error=mut_ddG['STD']
left=range(0,size(mut_ddG))
fig=plt.bar(left, ener, width=0.8, color="#ba00ff")
title= 'Alanine Scan'+title_1
plt.title(title)
plt.ylabel(r'ddG kcal mol$^{-1}$')
plt.xlabel('Residue')
# get axes handle
AXe = plt.gca()
# To get the current ticks
#mthd1
#xtic=AXe.get_xticks()
#or mthd2
#xtic=plt.getp(AXe,'xticks')

#==========
#   Text
#==========
#====== List RESN
# get xlim and ylim. add text after
xmax=plt.getp(plt.gca(),'xlim')[1]
ymax=plt.getp(plt.gca(),'ylim')[1]
yrange=plt.getp(plt.gca(),'ylim')[1]-plt.getp(plt.gca(),'ylim')[0]
for name in mut_ddG['resn']:
	plt.text(xmax*1.05,ymax,name,horizontalalignment='right',verticalalignment='top',fontsize=5)
	ymax=ymax-yrange/len(mut_ddG)

#====Label on top of graph
for i in range(0,len(mut_ddG)):
		if abs(ener[i])>0.4:
				plt.text(left[i]-xmax*0.005,ener[i]+ymax*0.05,mut_ddG['resn'][i],fontsize=7, bbox=dict(facecolor='orange', alpha=0.1))

#====== To set new ticks (pick every2)
n=2 # xticks interval
strt=2 # strting residue
ARR=array(left)
ARR=ARR[strt:len(ARR):n]
Name=mut_ddG['resn'][strt:len(mut_ddG):n]
plt.xticks(ARR,Name,size=7)
fig=plt.gcf()
fig.autofmt_xdate(ha='left',rotation=90)
plt.grid(True)
#plt.show()
plt.savefig('mut_ddG_bar.pdf', format='pdf')
##}}}
