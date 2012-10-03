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
      JOB="qsub 03_batch.sh"
      ;;
	m)
      JOB="qsub 03_mpi.sh"
      ;;
 	l)
      JOB="./03_local.sh &"
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

#############################
##### Sending the job
#############################

rm _MMPBSA_*
# Will not work on mac's bash thanks to different sed version
i=$(cd ../.. ; basename `pwd`)
sed -i -e "s|_NAME_|${i}|g" -e "s|_DIR_|`pwd`|g" 03_batch.sh
sed -i -e "s|_NAME_|${i}|g" -e "s|_DIR_|`pwd`|g" 03_mpi.sh
sed -i -e "s|_NAME_|${i}|g" -e "s|_DIR_|`pwd`|g" 03_local.sh
$JOB
