#!/bin/sh

# The purpose of this script is to organize all of the energy vectors 
# from an MMPBSA run into two files, one that has all of the bonded
# terms (bond, angle, dihedral, 1-4 electrostatics, 1-4 vdw) and the
# other that has all of the non-bonded terms (egb/epb, electrostatics
# and vdw) in another file.

# NOTE:
# If using a single trajectory, differences between bond lengths, angles,
# or dihedrals between the complex and receptor/ligand structures should
# cancel completely. This includes the BOND, ANGLE, DIHED, and 1-4 interactions.
# The calculation of delta G for single trajectory thus does not consider them 
# unless specified.


#######################################################################
#                                                                     #
#                User-defined variables. change as needed             #
#                                                                     #
#######################################################################

# _MMPBSA_ or _MMPBSA_mutant_ ?
prefix="_MMPBSA_"

# gb or pb mdouts? (lowercase only!)
solvent="gb"

# non-bonded output file
nbout="${solvent}_nonbond.dat"

# bonded output file (defaults to /dev/null, but you can put a name here)
bout="${solvent}_bond.dat"

#######################################################################
#                                                                     #
#                   Begin the actual script                           #
#                                                                     #
#######################################################################
# check if mdout.0 exist
exi=$(ls *.mdout.0 | wc -w)
if [ $exi -eq 0 ]; then
	#if mdout.0 does not exist
	ext=.mdout
else
	ext=.mdout.0
fi

if [ "`whoami`" = "root" ]; then
   echo "WARNING: DO NOT RUN THIS SCRIPT AS ROOT!"
   exit 1
fi

if [ ! -f ${prefix}complex_${solvent}${ext} ]; then
   echo "Error: Can't find ${prefix}complex_${solvent}${ext} ! MMPBSA"
   echo "       files are missing! Edit energyvectors.sh to specify"
   echo "       solvent model and/or prefix (_MMPBSA_ or _MMPBSA_mutant_)"
   exit 1
fi

# set caps-ed solvent
if [ "$solvent" = "gb" ]; then
   sol="GB"
else
   sol="PB"
fi

# remove old output files
echo "All files are listed in 3 columns: Complex  Receptor  Ligand"
echo "Removing ${nbout} and ${bout}..."
rm -f $nbout $bout .*evtmp 2>/dev/null
for i in Bond Angle Dihedral VDWAALS EEL E${sol} OnetoFourvdw ESURF;do
	rm ${i}.dat
done
echo "Files removed. Creating new files..."

# start with complex: store bond, angle, dihedral, 1-4 eel (eelof), 1-4 vdw
# (vdwof) vdw, eel, and esurf in .___.evtmp files

awk '$1=="BOND"{print $3}' ${prefix}complex_${solvent}${ext} > .bondc.evtmp
awk '$1=="BOND"{print $6}' ${prefix}complex_${solvent}${ext} > .anglec.evtmp
awk '$1=="BOND"{print $9}' ${prefix}complex_${solvent}${ext} > .dihedc.evtmp
awk '$1=="VDWAALS"{print $3}' ${prefix}complex_${solvent}${ext} > .vdwc.evtmp
awk '$1=="VDWAALS"{print $6}' ${prefix}complex_${solvent}${ext} > .eelc.evtmp
awk '$1=="VDWAALS"{print $9}' ${prefix}complex_${solvent}${ext} > .egbc.evtmp
awk '$1=="1-4"{print $4}' ${prefix}complex_${solvent}${ext} > .vdwofc.evtmp
awk '$1=="1-4"{print $8}' ${prefix}complex_${solvent}${ext} > .eelofc.evtmp
awk '$1=="1-4"{getline;print $3}' ${prefix}complex_${solvent}${ext} > .esurfc.evtmp

awk '$1=="BOND"{print $3}' ${prefix}receptor_${solvent}${ext} > .bondr.evtmp
awk '$1=="BOND"{print $6}' ${prefix}receptor_${solvent}${ext} > .angler.evtmp
awk '$1=="BOND"{print $9}' ${prefix}receptor_${solvent}${ext} > .dihedr.evtmp
awk '$1=="VDWAALS"{print $3}' ${prefix}receptor_${solvent}${ext} > .vdwr.evtmp
awk '$1=="VDWAALS"{print $6}' ${prefix}receptor_${solvent}${ext} > .eelr.evtmp
awk '$1=="VDWAALS"{print $9}' ${prefix}receptor_${solvent}${ext} > .egbr.evtmp
awk '$1=="1-4"{print $4}' ${prefix}receptor_${solvent}${ext} > .vdwofr.evtmp
awk '$1=="1-4"{print $8}' ${prefix}receptor_${solvent}${ext} > .eelofr.evtmp
awk '$1=="1-4"{getline;print $3}' ${prefix}receptor_${solvent}${ext} > .esurfr.evtmp

awk '$1=="BOND"{print $3}' ${prefix}ligand_${solvent}${ext} > .bondl.evtmp
awk '$1=="BOND"{print $6}' ${prefix}ligand_${solvent}${ext} > .anglel.evtmp
awk '$1=="BOND"{print $9}' ${prefix}ligand_${solvent}${ext} > .dihedl.evtmp
awk '$1=="VDWAALS"{print $3}' ${prefix}ligand_${solvent}${ext} > .vdwl.evtmp
awk '$1=="VDWAALS"{print $6}' ${prefix}ligand_${solvent}${ext} > .eell.evtmp
awk '$1=="VDWAALS"{print $9}' ${prefix}ligand_${solvent}${ext} > .egbl.evtmp
awk '$1=="1-4"{print $4}' ${prefix}ligand_${solvent}${ext} > .vdwofl.evtmp
awk '$1=="1-4"{print $8}' ${prefix}ligand_${solvent}${ext} > .eelofl.evtmp
awk '$1=="1-4"{getline;print $3}' ${prefix}ligand_${solvent}${ext} > .esurfl.evtmp

paste -d"   " .bondc.evtmp .bondr.evtmp > .evtmp
paste -d"   " .evtmp .bondl.evtmp > .bondt.evtmp
paste -d"   " .anglec.evtmp .angler.evtmp > .evtmp
paste -d"   " .evtmp .anglel.evtmp > .anglet.evtmp
paste -d"   " .dihedc.evtmp .dihedr.evtmp > .evtmp
paste -d"   " .evtmp .dihedl.evtmp > .dihedt.evtmp
paste -d"   " .vdwc.evtmp .vdwr.evtmp > .evtmp
paste -d"   " .evtmp .vdwl.evtmp > .vdwt.evtmp
paste -d"   " .eelc.evtmp .eelr.evtmp > .evtmp
paste -d"   " .evtmp .eell.evtmp > .eelt.evtmp
paste -d"   " .egbc.evtmp .egbr.evtmp > .evtmp
paste -d"   " .evtmp .egbl.evtmp > .egbt.evtmp
paste -d"   " .vdwofc.evtmp .vdwofr.evtmp > .evtmp
paste -d"   " .evtmp .vdwofl.evtmp > .vdwoft.evtmp
paste -d"   " .eelofc.evtmp .eelofr.evtmp > .evtmp
paste -d"   " .evtmp .eelofl.evtmp > .eeloft.evtmp
paste -d"   " .esurfc.evtmp .esurfr.evtmp > .evtmp
paste -d"   " .evtmp .esurfl.evtmp > .esurft.evtmp

# Time to compile everything:
echo "Complex Receptor Ligand" > Bond.dat
cat .bondt.evtmp >> Bond.dat
echo "Complex Receptor Ligand" > Angle.dat
cat .anglet.evtmp >> Angle.dat
echo "Complex Receptor Ligand" > Dihedral.dat
cat .dihedt.evtmp >> Dihedral.dat
echo "Complex Receptor Ligand" >> OnetoFourvdw.dat
cat .eeloft.evtmp >> OnetoFourvdw.dat
echo "Complex Receptor Ligand" > VDWAALS.dat
cat .vdwt.evtmp >> VDWAALS.dat
echo "Complex Receptor Ligand" > EEL.dat
cat .eelt.evtmp >> EEL.dat
echo "Complex Receptor Ligand" > E${sol}.dat
cat .egbt.evtmp >> E${sol}.dat
echo "Complex Receptor Ligand" >> ESURF.dat
cat .esurft.evtmp >> ESURF.dat

rm -f .*evtmp

# Creating a python script for computing delta and for plotting
python2.7 Calc_Delta.py ${sol}

# Create R script for computing delta
rm Calc_Delta.R; touch Calc_Delta.R
for i in Bond Angle Dihedral VDWAALS EEL E${sol} OnetoFourvdw ESURF;do
	echo "dat_${i} <- read.table(\"${i}.dat\",header=TRUE)" >> Calc_Delta.R
	echo "${i}_D <- dat_${i}[,1] - dat_${i}[,2] - dat_${i}[,3]" >> Calc_Delta.R
done
echo "SUM <- array(0,c(length(Bond_D),1))" >> Calc_Delta.R

# If separate trajectories are used, replace following line with commented one:
for i in VDWAALS EEL E${sol} ESURF;do
#for i in Bond Angle Dihedral VDWAALS EEL E${sol} OnetoFourvdw ESURF;do
	echo "SUM <- SUM + ${i}_D " >> Calc_Delta.R
	echo "mean(${i}_D) " >> Calc_Delta.R
done
echo "write.table(SUM,file = \"D_GBTOT\",quote=FALSE)" >> Calc_Delta.R

# Run R script
#R --vanilla < Calc_Delta.R
	

#echo "Done compiling energy vectors! non-bonded energies are in ${nbout}"
#echo "  bonded energies are in ${bout}!"
