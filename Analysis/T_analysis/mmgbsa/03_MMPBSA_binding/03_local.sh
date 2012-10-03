#!/bin/bash
#PBS -N _NAME_
#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -l walltime=50:00:00

############  DEFINE WORKING DIRECTORY   ####################
workdir=_DIR_
refdir="../../.."
#############################################################

#-------------------- Defining the nodes and cpus ----------#
#----------------------------- start the job -----------------------#

export AMBERHOME=/cluster/apps/x86_64/packages/amber11-patch-19/src

cd $workdir

STATUS_FILE=$workdir/STATUS
echo "START `date`" >> $STATUS_FILE
num_mdcrd=$(ls -1 $refdir| grep -P "prod\d+\.mdcrd\.gz" | wc -l)
num_mdcrd=$(($num_mdcrd*5000))
sed -i "s|_FRAME_|$num_mdcrd|g" mmpbsa.in
$AMBERHOME/bin/MMPBSA -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat -cp "$refdir/prot_0wat.parmtop" -rp "$refdir/rec.parmtop" -lp "$refdir/lig.parmtop" -y "$refdir/prod_vac.nc" > mmpbsa.out
sed -i "s|$num_mdcrd|_FRAME_|g" mmpbsa.in

echo "END `date`" >> $STATUS_FILE
