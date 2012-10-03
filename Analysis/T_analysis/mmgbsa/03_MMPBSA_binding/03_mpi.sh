#!/bin/bash
#PBS -N _NAME_
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -q ibqueue2
####PBS -l nodes=reindeer21
#PBS -l walltime=50:00:00

############  DEFINE WORKING DIRECTORY   ####################
workdir=_DIR_
refdir="../../.."
#############################################################

#-------------------- Defining the nodes and cpus ----------#
if [ ! -f $HOME/.mpd.conf ] ; then
       echo "FATAL! No $HOME/.mpd.conf"
       exit 1
     fi
     NODES=1   #<---- modify number of nodes here
     TOTAL_CPUS=$((NODES * 8))
     MPI_DIR=/cluster/apps/x86_64/packages/mvapich2-1.2p1-IB

    # Parse $PBS_NODEFILE, and generate $MACHINE_FILE
    MACHINE_FILE=$HOME/mpi/machine_file.$$
    sort < $PBS_NODEFILE | uniq > $MACHINE_FILE

    # Startup MPD, do a quick test before we proceed.
    $MPI_DIR/bin/mpdboot -n $NODES -f $MACHINE_FILE
    if [ $? -ne 0 ] ; then
      exit 1
    fi
    RESULT=`$MPI_DIR/bin/mpdtrace | wc -l`
    if [ "$RESULT" != "$NODES" ] ; then
      $MPI_DIR/bin/mpdallexit
      exit 1
    fi
#----------------------------- start the job -----------------------#
#start the job, shutdown MPD when complete.
#$MPI_DIR/bin/mpiexec -np $TOTAL_CPUS /cluster/apps/x86_64/bin/parallel-hello  <----- run your mpi program here 
#$MPI_DIR/bin/mpdallexit
rm -f $MACHINE_FILE

export AMBERHOME=/cluster/apps/x86_64/packages/amber12
cd $workdir

STATUS_FILE=$workdir/STATUS
echo "START `date`" >> $STATUS_FILE
num_mdcrd=$(ls -1 $refdir| grep -P "prod\d+\.mdcrd.*" | wc -l)
# Checks for folder accuracy
if [ "$num_mdcrd" -eq 0 ]; then
	echo "No mdcrd files found. Unable to determine number of frames."
	exit 1
fi
num_mdcrd=$(($num_mdcrd*5000)) #Total number of frames
sed -i "s|_FRAME_|$num_mdcrd|g" mmpbsa.in
$AMBERHOME/bin/MMPBSA.py.MPI -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat -cp "$refdir/prot_0wat.parmtop" -rp "$refdir/rec.parmtop" -lp "$refdir/lig.parmtop" -y "$refdir/prod_vac.nc" > mmpbsa.out
sed -i "s|$num_mdcrd|_FRAME_|g" mmpbsa.in
echo "END `date`" >> $STATUS_FILE
