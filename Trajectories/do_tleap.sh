#!/bin/bash

# greps away the water from prot.pdb first
# input: need prot.pdb in the folder
# output: prot_vac top/crd; prot_0wat top/crd; prot top/crd (has water box)
# if original prot.pdb has no waters, then the 0wat and vacs are the same

# list of output files:
# (0wat means no waters; vac means whatever the prot.pdb has, but with
# added hydrogens; nothing means water has been added)

#prot_0wat.pdb
#prot_vac.pdb
#prot_wat.pdb
#prot_0wat.inpcrd
#prot_0wat.parmtop
#prot_vac.inpcrd
#prot_vac.parmtop
#prot.inpcrd
#prot.parmtop


#export AMBERHOME=/cluster/apps/x86_64/packages/amber11-patch/src
#export AMBERHOME=/home/liuy/amber10

# override order temporarily to use the above amberhome's tleap instead of
# whatever is default
#export PATH=$AMBERHOME/bin:$PATH


# remove waters from pdb manually (tleap stupidly only removes 1 wat at a time)
# first few conditions retain only non-HOH and non-WAT lines, while keeping TER
# and END cards. uniq removes consecutive TER cards (any atom should not be
# identical from line to line).
awk '((substr($0,18,3) != "HOH") && (substr($0,18,3) != "WAT")) || ($1 == "TER") || ($1 == "END") ' prot.pdb | uniq > prot_0wat.pdb

tleap -f leap.in > myleap.log

