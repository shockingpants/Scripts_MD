#!/opt/local/bin/python2.7
# fasta_parse.py
from strtup import *
import sys
from Bio import SeqIO
import numpy
import matplotlib.pyplot as plt

cd('/Volumes/HDD/teojy/simulations/MDMX_P53/ori_file')
target=open('prot_seq.fasta','r')
j=0
name={}
sequence={}
for i in SeqIO.parse(target, "fasta"):  # It is a list of classes, stored from the generator
    j=j+1
    name[j] = i.id
    sequence[name[j]]=i._seq._data

## Set offset between actual file index and aminoacid index
offset={}
offset['human_mdm4']=23
offset['mouse_mdm4']=23
offset['zebrafish_mdm4']=23
offset['human_mdm2']=24

pathdir=['human_mdmx_p53', 'human_mdmx_p53_2', 'human_mdmx_p53_3', 'mouse_mdmx_p53', 'mouse_mdmx_p53_2', 'mouse_mdmx_p53_3', 'zebrafish_mdmx_p53', 'zebrafish_mdmx_p53_2', 'zebrafish_mdmx_p53_3']

for k in range(0,9):
    cd('/Volumes/HDD/teojy/simulations/MDMX_P53/'+pathdir[k])
    mkdir('ala_scan')
    cd('ala_scan')
    numb = (k+3)/3
    for aa in range(0,len(sequence[name[numb]])):
        dir_name=str(aa+offset[name[numb]])+'_'+sequence[name[numb]][aa]
        mkdir(dir_name)
    
    
    
    
