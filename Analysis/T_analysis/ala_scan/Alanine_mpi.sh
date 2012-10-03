#!/bin/bash
#PBS -N _NAME_
#PBS -l nodes=1:ppn=8
#PBS -j oe
#PBS -q ibqueue2
####PBS -l nodes=reindeer21
#PBS -l walltime=80:00:00

############  DEFINE WORKING DIRECTORY   ####################
workdir=_PWD_
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

$AMBERHOME/bin/MMPBSA.py.MPI -O -i mmpbsa_mpi.in -o FINAL_RESULTS_MMPBSA.dat -cp ../../prot_0wat.parmtop -rp ../../rec.parmtop -lp ../../lig.parmtop -y ../../prod_vac.nc -mc mut.parmtop -mr mut_rec.parmtop > mmpbsa.out 
                         
echo "END `date`" >> $STATUS_FILE

# remove trajectory
rm *.mdcrd*


