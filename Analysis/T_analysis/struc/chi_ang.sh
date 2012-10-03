#!/bin/bash
# chi_ang.sh
# Extract chi angles
# Called from chi_angles.py

################################
####	Option Parsing
################################
##{{{
#usage()
###{{{
#{
#cat << EOF
#Usage: $0 [options]
#EG   : $0 -i -p -y
#
#OPTIONS:
#	-h	show this help message and exit
#	-i	input scripts
#	-p	prmtop file
#	-d	pdb file
#	-o	output file
#	-y	trajectory file
#EOF
#}
###}}}
#
#if [ -z $1 ]; then
#	echo "Jon: Please choose an option/argument!"
#	usage
#	exit 1
#fi
#INPUT=
#OUTPUT=
#PRMTOP=
#PDB=
#TRAJ=
#while getopts "h:i:p:d:o:y" opt; do
#  case $opt in
#    h)
#      usage
#	  exit 1
#	  ;;
#	i)
#      INPUT=$OPTARG
#      ;;
# 	p)
#      PRMTOP=$OPTARG
#	  ;;
#	d)
#      PDB=$OPTARG
#	  ;;
#	o)
#      OUTPUT=$OPTARG
#      ;;
#	y)
#      TRAJ=$OPTARG
#      ;;
#
#    \?) # This is returned when the option does not correspond to hipoy
#      echo "Invalid option: -$OPTARG" >&2
#      exit 1
#      ;;
#    :) # This is returned when an argument is NOT specified for an option
#      echo "Option -$OPTARG requires an argument." >&2
#      exit 1
#      ;;
#  esac
#done
#
#if [[ -z $PARMTOP ]] || [[ -z $TRAJ ]]; then
#     usage
#     exit 1
#fi
#
##}}}

################################
####       Start up
################################
echo "Loading pdb into vmd...      "   `date`
wdir=`pwd`
# cd to right folder
cd ../../
echo `pwd`
/Applications/VMD_1.9.1.app/Contents/MacOS/startup.command
cd $wdir

################################
####     Generate Chi Angles
################################
echo "Generating Chi Angles...      "   `date`
mkdir $FOLD
cd $FOLD
echo `pwd`
g_chi -s ../../../prot_0wat.pdb -f ../../../prod_vac.trr -maxchi 1 -all


################################
####	Extract Chi Angles
################################
echo "Extracting Chi Angles...      "  `date`
ls -1
rm chi_index.txt
for i in $(ls -1 chi*|grep "chi1[A-Z]\{3\}.*" | gsed "s/\(chi1[A-Z]\{3\}\)\([0-9]\+\).xvg/\2/g"|sort -n); do
	echo $i
	filename=$(ls -1 chi*| grep -P "chi1[A-Z]{3}$i\..*")
	echo $i $filename >> chi_index.txt
done


