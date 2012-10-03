#!/bin/bash
# you can specify your own range
# choose ibq or batch or local
# Options
# -b batch
# -m mpi
# -l local

###############################
##### User Input
###############################

if [ -z $1 ]; then
		echo "Jon: Please choose an option/argument!"
		exit 1
fi

while getopts ":bml" opt; do
  case $opt in
    b)
      JOB="qsub Alanine_batch.sh"
      ;;
	m)
      JOB="qsub Alanine_mpi.sh"
      ;;
 	l)
      JOB="./Alanine_local.sh &"
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

##========= Change the for loop below accordingly
#RANGE=(99_Y)


#############################
##### Sending the job
#############################

rm _MMPBSA_*

# This command finds all folders except self. Can replace with RANGE
for i in $(ls -1| grep '^[0-9][0-9]*_[A-Z]'); do
#for i in $(find . -type d -not -name .); do
#for i in ${RANGE[@]}; do
	i=$(echo $i| sed "s|^\.\/||")
    cd ${i}
    # the `pwd` only works for certain sed.
	cp ../Alanine_batch.sh .
	cp ../Alanine_mpi.sh .
	cp ../Alanine_local.sh .
	#cp ../nc_conv.sh .
	cp ../mmpbsa.in .
	cp ../mmpbsa_mpi.in .
    sed -i -e "s|_NAME_|${i}|g" -e "s|_PWD_|`pwd`|g" Alanine_batch.sh
    sed -i -e "s|_NAME_|${i}|g" -e "s|_PWD_|`pwd`|g" Alanine_mpi.sh
	sed -i -e "s|_NAME_|${i}|g" -e "s|_PWD_|`pwd`|g" Alanine_local.sh
	chmod 777 Alanine_batch.sh
    chmod 777 Alanine_mpi.sh
	chmod 777 Alanine_local.sh
	#chmod 777 nc_conv.sh
	$JOB
	cd ..
done
jobs -l
