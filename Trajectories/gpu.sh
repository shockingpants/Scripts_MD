#!/bin/bash
#PBS -j oe
#PBS -q cuda 
#PBS -N _NAME_


###### preamble starts #######
# determine GPU number by doing a simple mem_test 
NUM_DEVICE=`/usr/bin/nvidia-smi -L | wc -l`


for I in `seq 0 $((NUM_DEVICE - 1))` 
do
/cluster/apps/x86_64/packages/utils/cuda_memtest --stress --num_passes 1 --num_iterations 1 --device $I > /dev/null 2>/dev/null
if [ "$?" -eq "0" ]; then
GPU_NUM=$I
break
fi
done
# $GPU_NUM is now "--device $I", dont change value of $GPU_NUM as inputs to cuda-fied binaries unless you know what you are doing !!!
##### preamble ends ########


sleep 1

#################################
# User defined part starts here #
#################################

############  DEFINE WORKING DIRECTORY   ####################
WORKDIR=_DIR_
STATUS=$WORKDIR/STATUS
#############################################################

##########################
##### Functions
##########################
# compress mdcrd function
compr(){
echo $1
if [ -f $1 -a ! -f $1.gz ] 
then
echo $1  `date`>> compress.log
echo "Compressing $1..."
gzip 9 $1
elif [ ! -f $1.gz ] 
then
echo "$1 is already compressed."  `date`>> compress.log
else
echo $1 does not exist. `date`>> compress.log
fi
}

##########################
####  Actual Code
##########################

export AMBERHOME=/cluster/apps/x86_64/packages/amber11-patch-19/
export PATH=$AMBERHOME/bin:$PATH:/cluster/apps/x86_64/packages/cuda/bin
RUN_CUDA="$AMBERHOME/bin/pmemd.cuda -gpu $GPU_NUM -O"

cd $WORKDIR
echo "START `date`" >> $STATUS


$RUN_CUDA -i min1.in  -o min1.out  -p prot.parmtop -c prot.inpcrd -r min1.rst  -ref prot.inpcrd </dev/null
echo -e "END_MD_MIN1\t`date`" >> $STATUS
$RUN_CUDA -i min2.in  -o min2.out  -p prot.parmtop -c min1.rst    -r min2.rst  -ref min1.rst </dev/null
echo -e "END_MD_MIN2\t`date`" >> $STATUS
$RUN_CUDA -i min3.in  -o min3.out  -p prot.parmtop -c min2.rst    -r min3.rst </dev/null
echo -e "END_MD_MIN3\t`date`" >> $STATUS
$RUN_CUDA -i heat1.in -o heat1.out -p prot.parmtop -c min3.rst    -r heat1.rst -x heat1.mdcrd </dev/null
echo -e "END_MD_HEAT1\t`date`" >> $STATUS
$RUN_CUDA -i heat2.in -o heat2.out -p prot.parmtop -c heat1.rst   -r heat2.rst -x heat2.mdcrd </dev/null
echo -e "END_MD_HEAT2\t`date`" >> $STATUS
$RUN_CUDA -i heat3.in -o heat3.out -p prot.parmtop -c heat2.rst   -r heat3.rst -x heat3.mdcrd </dev/null
echo -e "END_MD_HEAT3\t`date`" >> $STATUS
$RUN_CUDA -i prod.in -o prod1.out -p prot.parmtop -c heat3.rst   -r prod1.rst -x prod1.mdcrd </dev/null
echo -e "END_MD1\t`date`" >> $STATUS
compr prod1.mdcrd &
$RUN_CUDA -i prod.in -o prod2.out -p prot.parmtop -c prod1.rst   -r prod2.rst -x prod2.mdcrd </dev/null
echo -e "END_MD2\t`date`" >> $STATUS
compr prod2.mdcrd &
$RUN_CUDA -i prod.in -o prod3.out -p prot.parmtop -c prod2.rst   -r prod3.rst -x prod3.mdcrd </dev/null
echo -e "END_MD3\t`date`" >> $STATUS
compr prod3.mdcrd &
$RUN_CUDA -i prod.in -o prod4.out -p prot.parmtop -c prod3.rst   -r prod4.rst -x prod4.mdcrd </dev/null
echo -e "END_MD4\t`date`" >> $STATUS
compr prod4.mdcrd &
$RUN_CUDA -i prod.in -o prod5.out -p prot.parmtop -c prod4.rst   -r prod5.rst -x prod5.mdcrd </dev/null
echo -e "END_MD5\t`date`" >> $STATUS
compr prod5.mdcrd &
$RUN_CUDA -i prod.in -o prod6.out -p prot.parmtop -c prod5.rst   -r prod6.rst -x prod6.mdcrd </dev/null
echo -e "END_MD6\t`date`" >> $STATUS
compr prod6.mdcrd &
$RUN_CUDA -i prod.in -o prod7.out -p prot.parmtop -c prod6.rst   -r prod7.rst -x prod7.mdcrd </dev/null
echo -e "END_MD7\t`date`" >> $STATUS
compr prod7.mdcrd &
$RUN_CUDA -i prod.in -o prod8.out -p prot.parmtop -c prod7.rst   -r prod8.rst -x prod8.mdcrd </dev/null
echo -e "END_MD8\t`date`" >> $STATUS
compr prod8.mdcrd &
$RUN_CUDA -i prod.in -o prod9.out -p prot.parmtop -c prod8.rst   -r prod9.rst -x prod9.mdcrd </dev/null
echo -e "END_MD9\t`date`" >> $STATUS
compr prod9.mdcrd &
$RUN_CUDA -i prod.in -o prod10.out -p prot.parmtop -c prod9.rst   -r prod10.rst -x prod10.mdcrd </dev/null
echo -e "END_MD10\t`date`" >> $STATUS
compr prod10.mdcrd &
#$RUN_CUDA -i prod.in -o prod11.out -p prot.parmtop -c prod10.rst   -r prod11.rst -x prod11.mdcrd </dev/null
#echo -e "END_MD6\t`date`" >> $STATUS
#$RUN_CUDA -i prod.in -o prod12.out -p prot.parmtop -c prod11.rst   -r prod12.rst -x prod12.mdcrd </dev/null
#echo -e "END_MD7\t`date`" >> $STATUS
#$RUN_CUDA -i prod.in -o prod13.out -p prot.parmtop -c prod12.rst   -r prod13.rst -x prod13.mdcrd </dev/null
#echo -e "END_MD8\t`date`" >> $STATUS
#$RUN_CUDA -i prod.in -o prod14.out -p prot.parmtop -c prod13.rst   -r prod14.rst -x prod14.mdcrd </dev/null
#echo -e "END_MD9\t`date`" >> $STATUS
#$RUN_CUDA -i prod.in -o prod15.out -p prot.parmtop -c prod14.rst   -r prod15.rst -x prod15.mdcrd </dev/null
#echo -e "END_MD10\t`date`" >> $STATUS
