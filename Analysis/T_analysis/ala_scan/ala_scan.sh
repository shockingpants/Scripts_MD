#!/bin/bash
# Input prot_0wat.pdb
# Start in ala_scan folder 
# Output Folders with pdb and parmtop of mut_{com, lig, rec}


#####################
#### Creating fasta files and setting parameters
#####################

read -p "Have you changed OFFSET, FIRST and END (y/n): " ans
if [[ ! $ans == "y" ]]; then
	echo Please enter a right option.
	echo Please change the variables.
	exit
fi

cd ..
PDB_file='prot_0wat'
OFFSET=_offset_
pdb2fasta "$PDB_file.pdb" $OFFSET > "$PDB_file.fasta"

#============== Actual resi in literature (aka index + offset)
# Choose one, uncomment where necessary
# EITHER
FIRST=_start_
END=_end_
RANGE=($(seq $FIRST $END))
# OR
#RANGE="manual input in the form of an array"
#RANGE

SIZE=${#RANGE[@]}



#########################
#### Parsing fasta files
#########################

#===== Get sequence (This combines several lines into 1)
sequ=$(egrep -v ">" "$PDB_file.fasta")
sequ=($(echo $sequ | sed "s| ||g"))


#===== Get number of aa
num_aa=$(echo $sequ | wc -c)
num_aa=$((num_aa-1))

###===== Get offset
##offset=$(grep '>' lig.fasta | awk '{print $6}')


########################
#### Check for other files
########################
for i in prod_vac.nc lig.pdb lig.parmtop lig.inpcrd rec.pdb rec.parmtop rec.inpcrd
do
	if [ ! -f $i ]
	then
		echo "$i does not exist"
		exit
	fi
	echo
done


########################
##### Making Directories, transfer files and creating parmtop
########################
cd ala_scan

for i in $(seq 0 $((SIZE-1))); do
    #===================
    # Residue check
    #===================
    # index on sequence
    IND=$((${RANGE[$i]}-$OFFSET))
    # Use index to call the residue name
    RESN=${sequ[0]:$((IND-1)):1}
    if [[ $RESN == "X" ]] || [[ $RESN == "P" ]] || [[ $RESN == "G" ]] || [[ $RESN == "A" ]]; then
        continue
    fi

    
    #==================
    # Make Directory and text file 
    #==================
	echo "${RANGE[$i]}_$RESN"
    dir_name=$(echo "${RANGE[$i]}_$RESN")
    mkdir -p $dir_name
	echo $dir_name >> folder_list.txt
    cd $dir_name
    
    #================
    # Creating the mutant 
    #================
    cp ../../prot_0wat.pdb .
    # lines to delete
    awk ' $3 != "C" && $3 != "C" && $3 != "N" && $3 != "CA" && $5 == '$IND' {print $1"_"$2" "}' prot_0wat.pdb > del.pdb
    cp prot_0wat.pdb mut.pdb
    OMIT=(`cat del.pdb`)
    rm del.pdb
    for i in ${OMIT[@]}; do
        i=$(echo $i | sed 's|_| *|')
        egrep -v "$i " mut.pdb > temp.pdb
        mv temp.pdb mut.pdb
    done
    sed -i '' -e 's/ [[:alpha:]][[:alpha:]][[:alpha:]]\( *'$IND' \)/ ALA\1/' mut.pdb
    
    #==============
    # Split files into different PDB
    #==============
    $SCRIPTS/splitpdb.py mut_rec.pdb mut_lig.pdb nil.pdb < mut.pdb
    rm nil.pdb
    
    #==============
    # Create Parmtops and stuff
    #==============
    tleap -f ../leap_mutreclig.in > leap_mutreclig.log
    
    #==============
    # Copy essential files
    #==============
    cp ../mmpbsa.in .
    cp ../Alanine_local.sh .
	cp ../Alanine_mpi.sh .
	cp ../Alanine_batch.sh .
	# cp ../nc_conv.sh .
        
    cd ..
done


#======== sort folder_list
cat folder_list.txt | sort -n > temp.txt
mv temp.txt folder_list.txt




