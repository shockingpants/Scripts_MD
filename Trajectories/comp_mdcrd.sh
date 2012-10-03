#!/bin/bash
#PBS -N _NAME_
#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -l walltime=50:00:00

# Run this only if batch is free!

############  DEFINE WORKING DIRECTORY   ####################
workdir=_DIR_
#############################################################

#-------------------- Defining the nodes and cpus ----------#
#----------------------------- start the job -----------------------#
cd $workdir

#Compression 9
for j in $(ls -1|grep '.*\.mdcrd$'); do
        # Checks for existence of .gz
        if [ -f ${j} -a ! -f ${j}.gz ]; then
                echo $j >> compress.log
                echo "Compressing ${j}"
                gzip 9 ${j}
        fi
done

