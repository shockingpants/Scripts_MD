#!/bin/bash
# This helps prepare files for use in MD for the MDMX_P53 project.
# Please switch on macfusion to allow sync
# Execute from within folder containing the desired PDB file.
###################################
####    Set up paths
###################################
##{{{
workdir=$PWD
scriptdir='/Volumes/HDD/teojy/Scripts' # Where the scripts are stored
#santadir='/Volumes/Santa' # Virtual disk of santa
##}}}
###################################
####     Set up directories
###################################
##{{{
# Check if script folder exist and copy files over

if [ -d "$scriptdir" ] # Use -f if checking for files
then
    cp -r "$scriptdir/Trajectories/" .  # -r for copying directories
    cp -r "$scriptdir/Analysis/T_Analysis" .
    cp "$scriptdir/splitpdb.py" .
    # Change permission
    chmod 700 do_tleap.sh
    chmod 700 gpu.sh
    chmod 700 ./mmgbsa/03*/03.sh
else
    echo "$scriptdir does not exist"
    exit
fi
##}}}
###################################
####    Setting up protein
###################################
# Remove H and OXT from pdb file
if [[ $(cat prot.pdb | grep ' H *$' | wc -l) -ne 0 ]]; then
cat prot.pdb | grep -v ' H *$' > temp.pdb
mv prot.pdb prot_H.pdb
mv temp.pdb prot.pdb
fi

if [[ $(cat prot.pdb | grep 'OXT' | wc -l) -ne 0 ]]; then
cat prot.pdb | grep -v 'OXT' > temp.pdb
mv prot_H.pdb prot_H_OXT.pdb
mv temp.pdb prot.pdb
fi


# Create a pdb file stripped of water
awk '((substr($0,18,3) != "HOH") && (substr($0,18,3) != "WAT")) || ($1 == "TER") || ($1 == "END") ' prot.pdb | uniq > prot_0wat.pdb

# Extracting fasta etc information after water is extracted
./resname.py


# Generate parmtop and inpcrd files
./do_tleap.sh

# Check for fatal errors
if [[ $(cat myleap.log | grep -i 'fatal' | wc -l) -ne 0 ]]; then
cat myleap.log | grep -i 'fatal' >> FATAL.log
echo THERE ARE FATAL ERRORS! CHECK myleap.log or FATAL.log
fi

# Edit min1.in file to reflect number of residues
# Assumes that there is a chain number
num_res=$(cat prot.pdb | awk ' $1=="ATOM" {print $4, $6}' | uniq | wc -l
)
num_res=$(($num_res+10))
sed -i '' -e "s|_RES_|$num_res|" min1.in
echo Check min1.in to see if res number is absurd!
echo Has to do with presence of chain letter or not


# Creates the receptor and ligand file
./Processes/split_pdb.py prot_0wat.pdb -o rec.pdb lig.pdb 

tleap -f leap_reclig.in >> myleap.log
date >> myleap.log
