#!/bin/bash
# Converts mdcrd files to an nc file.
# This will help save space.
# Remove heavy ions and water
# Options: 
# -a all .mdcrd files in the folder
# -s select individual .mdcrd file

while getopts ":a:s:" opt; do
		case $opt in
		a) 
			# Stores all existing .mdcrd files into frames

		;;

		s)
			# Stores the .mdcrd files into frames
			echo 'Please choose the .mdcrd files you would like to process into an nc trajectory file'
			ls *.mdcrd
			i=0
			j=0
			while [ $i -ne 1 ];
				do
   				j=$(($j+1))
    			read -p "Which frame would you like to extract?  " frame[$j]
					# Check if input is null
   					if [ -z ${frame[$j]} ]; then
        				i=1
        				break
   					fi
			done

		;;

		\?)
		 echo "Invalid option: -$OPTARG" >&2
      	exit 1
      	;;

		:)
	    echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
		esac
done



		
#echo 'Please choose the .mdcrd files you would like to process into an nc trajectory file'
#ls *.mdcrd
#i=$((0))
#j=$((0))
#while [ $i -ne $((1)) ];
#do
#    j=$(($j+1))
#    read -p "Which frame would you like to extract?  " frame[$j]
#    if [ -z ${frame[$j]} ]; then
#        i=1
#        break
#    fi
#done

END_RES=`awk "/Na|Cl|WAT/{exit};1" prot_wat.pdb | awk '/^ATOM/{print $5}' |     tail -n1`

ptraj prot.parmtop << _END_ > ptraj_stripwat_out2nc.out
#for k in $(seq 1 1 $j); do
#trajin ${frame[$k]}
#done
## if restarting from a previous simulation: (assume we are only keeping first 50 frames
## trajin prod1.mdcrd 1 50 1

trajin prod1.mdcrd
trajin prod2.mdcrd
trajin prod3.mdcrd
trajin prod4.mdcrd
trajin prod5.mdcrd

## This strips out the water and heavy ions
strip  :WAT
strip  :Cl-
strip  :Na+
reference prot_0wat.pdb

center :1-${END_RES}
image center

rms reference mass out ptraj_rms.dat @CA time 2.0
trajout prod_vac.nc netcdf nobox

go
_END_


# get a CA-only trj
ptraj prot_0wat.parmtop << _END_ > ptraj_CA_out2nc.out
trajin  prod_vac.nc
strip !@CA
trajout prod_CA.nc netcdf nobox
go
_END_

# spit out CA pdb
ptraj prot_0wat.parmtop << _END_ > ptraj_CA_out2pdb1.out
trajin  prod_vac.nc 1 1
strip !@CA
trajout prod_CA.pdb pdb
go
_END_

