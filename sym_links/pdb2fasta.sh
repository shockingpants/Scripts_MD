#!/bin/bash
#######################################
#############  Extract resn and resid
# http://cupnet.net/pdb2fasta/
# https://bitbucket.org/pierrepo/pdb2fasta/changeset/2a7191edac8f#chg-pdb2fasta.sh
#######################################

# We get the 3 letter sequence by selecting for first column == atom, 3rd column == alpha carbon and 5th column = desired chain
# To remove /n for each line (aka transpose) try ===> sed ':a;N;$!ba;s/\n/ /g' or ==> tr
# This does not take multiple chain into account. Read bit bucket link for more details

#============
# Usage
#============
# ./pdb2fasta.sh *.pdb offset > *.fasta

#=================
# Functions
#=================

# Extract residue from sequence lines
extract() {
    # $1 here refers to the first input following this function
    # -v inputs a variable. ch can then be called in the awk script. Its redundant here, but for reference here
    awk -v ch=$1 '/^ATOM/ && $3 == "CA" || $3 == "CH3"{print $4}'
}

# convert newline by space
# for sed lowers this works too: sed ':a;N;$!ba;s/\n/ /g'
remove_newline() {
    tr '\n' ' ' 
}


# convert 3-letter residue code to 1-letter
convert_aa() {
    sed 's/ALA/A/g;s/CYS/C/g;s/ASP/D/g;s/GLU/E/g;s/PHE/F/g;s/GLY/G/g;s/HIS/H/g;s/ILE/I/g;s/LYS/K/g;s/LEU/L/g;s/MET/M/g;s/ASN/N/g;s/PRO/P/g;s/GLN/Q/g;s/ARG/R/g;s/SER/S/g;s/THR/T/g;s/VAL/V/g;s/TRP/W/g;s/TYR/Y/g;s/HIE/H/g;s/ACE/X/g;s/NME/X/g'
}
# We get LYS LYS ALA VAL ILE ASN ...

# remove space between residues
remove_space() {
    sed 's/ //g' 
}

# split fasta sequence at 60 characters (easier to read than 80)
split_60() {
    fold -w 60
}

#======================
# Run the code
#======================
# First argument after the script. Typically prot_0wat.pdb
PDB_file=$1
# Second argument is the offset from the index
# Eg. What does the first residue correspond to
offset=$2
sequence=$(cat $PDB_file | extract | remove_newline | convert_aa | remove_space)
if [[ -n $sequence ]]; then
    size=$(echo $sequence | wc -c)
    size=$(($size-1))
    if [[ -n $offset ]]; then
        echo ">$PDB_file | $size aa | $offset offset"
    else
        echo ">$PDB_file | $size aa | 0 offset"
    fi
    echo $sequence | split_60
fi