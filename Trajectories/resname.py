#!/opt/local/bin/python2.7
# Extract filename and indexes into various formats
# fasta, genbank, txt
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
rpt='1_Repeat'
frst_rpt='../'+rpt+'/' #Folder containing the first repeat
trj_folder='./'

############################################
####     Classes and functions
############################################
##{{{
# 1 to 3
aa_codes = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU',
                    'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE',
                    'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
                    'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER',
                    'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR',
					'd':'HID', 'e':'HIE', 'p':'HIP', 'X':'XXX',}
# 3 to 1
inverse_aa_codes = dict([(three, one) for one,three in aa_codes.items()])

# from GromacsWrapper utilities
def conv_aa(x):
	"""Converts between 3-letter and 1-letter amino acid codes."""
	if len(x) == 1:
		if x not in aa_codes:
			return 'XXX'
		else:	
			return aa_codes[x]
	elif len(x) == 3:
		if x not in inverse_aa_codes:
			return 'X'
		else:
			return inverse_aa_codes[x]
	else:
		raise ValueError("Can only convert 1-letter or 3-letter amino acid codes,not %r" % x)
##}}}
############################################
####   Extracting original pdb file
############################################
##{{{
# Check if frst_rpt exist. If it does, copy the files over. Simple
import shutil
if os.path.basename(os.getcwd())!=rpt: 	# Copy file over
	if os.path.isfile(frst_rpt+'resname.txt'):   
		# Checks if file exist in frst_rpt. guard against older trajectories
		shutil.copy(frst_rpt+'resname.txt','.') # Text file with names and offset
		shutil.copy(frst_rpt+'resname.fasta','.') # Fasta file
		shutil.copy(frst_rpt+'resname.gbk','.') # Genbank file
	else:
		raise Exception('File does not exist. Please run resname.py in '+frst_rpt+'.')

	sys.exit('Copying over from '+frst_rpt)


print 'Extracting data from prot.pdb'
#=======================================
# 		Resname.txt
#=======================================
##{{{
#--------  R -----------
#b3d=importr('bio3d')
#pdb=b3d.read_pdb(trj_folder+"prot.pdb") # Crystal Structure
#ca_inds2=b3d.atom_select(pdb, "calpha").rx2('atom') # Crystal Structure with calpha only
#ca_seq=b3d.aa321(pdb.rx2('atom').rx(ca_inds2,'resid'))
#ca_resno=pdb.rx2('atom').rx(ca_inds2,'resno')
#ca_names=R.paste(ca_seq,ca_resno) # Y100 for example

#-------- MDA ----------
uni=MDA.Universe('prot_0wat.pdb')
atom=uni.selectAtoms('all')
ca_seq=atom.resids()
ca_names=atom.resnames()
ca_letter=map(conv_aa,ca_names)
# index starts from zero
with open('resname.txt','w') as f:
	f.write('Total # of residues = {0:d} \n'.format(len(ca_seq))+'\n')
	#f.write('Offset = ['+', '.join(map(str,np.array(ca_seq)-np.arange(len(ca_seq))))+']\n') # Offset is based on difference between actual resid and python index starting from zero
	f.write('\t'.join(['ind','resid','offset','three_aa','id_name1','id_name3','name1_id','name3_id'])+'\n')
	for ind,i in enumerate(ca_seq):
		f.write('\t'.join(map(str,[ind,i,ind-i,ca_names[ind],str(i)+ca_names[ind],str(i)+ca_letter[ind],ca_names[ind]+str(i),ca_letter[ind]+str(i)]))+'\n')
##}}}
#=======================================
# 		Resname.fasta
#=======================================
##{{{
# Creatin fasta file
n=60 # Number of letter per line
with open('resname.fasta','w') as f:
	f.write('> pdb|{0:s}'.format(uni.filename)+'\n')
	rows=len(ca_seq)/n+1
	for i in xrange(rows):
		if i != (rows-1):
			f.write(''.join(map(str,ca_letter[i*n:i+1*n]))+'\n')
		else:
			f.write(''.join(map(str,ca_letter[i*n:])))
##}}}
#=======================================
#		Resname.fasta2
#=======================================
##{{{
n=50 # Number of letters per line
nn=10# 10 per block
ca_space=[' ' if ind%nn != 0 else str(ind+1) for ind,i in enumerate(ca_letter)] # spaces and index
ca_list1=[]
ca_list2=[]
for i in range(len(ca_letter))[::nn]:
	ca_list1.append(''.join(ca_letter[i:i+nn]))
	ca_list2.append(''.join(ca_space[i:i+nn]))

with open('resname.fasta2','w') as f:
	f.write('> pdb|{0:s}'.format(uni.filename)+'\n')
	rows=len(ca_seq)
	for i in xrange(rows):
		if i != (rows-1):
			f.write('\t'.join(ca_list2[i*(n/nn):(i+1)*(n/nn)])+'\n')
			f.write('\t'.join(ca_list1[i*(n/nn):(i+1)*(n/nn)])+'\n\n')
		else:
			f.write('\t'.join(ca_list2[i*(n/nn):])+'\n')
			f.write('\t'.join(ca_list1[i*(n/nn):]))
##}}}
#=======================================
# 		Resname.gbk
#=======================================
##{{{
from Bio import SeqIO,Alphabet
input_handle = open("resname.fasta", "rU")
output_handle = open("resname.gbk", "w")
sequences = SeqIO.parse(input_handle, "fasta")
for seq in sequences:
	seq.seq.alphabet=Alphabet.generic_protein
	seq.name = seq.name.split('|')[1]
	SeqIO.write((seq,), output_handle, "genbank")

output_handle.close()
input_handle.close()
##}}}
##}}}
		
	
	

