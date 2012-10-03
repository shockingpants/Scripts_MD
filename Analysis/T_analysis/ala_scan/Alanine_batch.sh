#!/bin/bash
#PBS -N _NAME_
#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -l walltime=50:00:00

############  DEFINE WORKING DIRECTORY   ####################
workdir=_PWD_
#############################################################

#-------------------- Defining the nodes and cpus ----------#
#----------------------------- start the job -----------------------#

export AMBERHOME=/cluster/apps/x86_64/packages/amber11-patch-19/src

cd $workdir

STATUS_FILE=$workdir/STATUS
echo "START `date`" >> $STATUS_FILE

$AMBERHOME/bin/MMPBSA -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat -cp ../../prot_0wat.parmtop -rp ../../rec.parmtop -lp ../../lig.parmtop -y ../../prod_vac.nc -mc mut.parmtop -mr mut_rec.parmtop > mmpbsa.out 
                         
echo "END `date`" >> $STATUS_FILE

chmod 777 nc_conv.sh
./nc_conv.sh
