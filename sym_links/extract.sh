#!/bin/bash
#PBS -N _NAME_
#PBS -l nodes=1:ppn=1
#PBS -q batch
#PBS -l walltime=10:00:00

############  DEFINE WORKING DIRECTORY   ####################
workdir=_DIR_
#############################################################

#-------------------- Defining the nodes and cpus ----------#
#----------------------------- start the job -----------------------#


# cd $workdir
# For lines with Na, Cl or WAT, exit, if not, print (1" means print) 
END_RES=`awk "/Na|Cl|WAT/{exit};1" prot_wat.pdb | awk '/^ATOM/{print $5}' | tail -n1`

#/cluster/apps/x86_64/packages/amber11-patch/src/bin/ptraj prot.parmtop << _END_ > ptraj_stripwat_out2nc.out

#====== Create NETcdf files
# Finds all prod*.mdcrd(.gz) files
for i in $(ls -1| grep -P '^prod\d+\.mdcrd*'); do
	echo "trajin ${i}" >> TEM.in
done

# Sorts the mdcrd files
cat TEM.in | sort -td -k2 -n > TEM2.in # NOTE: this works only if name follows convention prod###
mv TEM2.in TEM.in

## if restarting from a previous simulation: (assume we are only keeping first 50 frames
## trajin prod1.mdcrd 1 50 1

echo "
strip  :WAT
strip  :Cl-
strip  :Na+
reference prot_0wat.pdb

center :1-${END_RES}
image center

rms reference mass out ptraj_rms.dat @CA time 2.0
trajout prod_vac.nc netcdf nobox

go
" >> TEM.in

ptraj prot.parmtop < TEM.in > ptraj_stripwat_out2nc.out
rm TEM.in

# get a CA-only trj
#/cluster/apps/x86_64/packages/amber11-patch/src/bin/ptraj prot_0wat.parmtop << _END_ > ptraj_CA_out2nc.out
ptraj prot_0wat.parmtop << _END_ > ptraj_CA_out2nc.out

trajin  prod_vac.nc
strip !@CA
trajout prod_CA.nc netcdf nobox
go
_END_

# spit out CA pdb
#/cluster/apps/x86_64/packages/amber11-patch/src/bin/ptraj prot_0wat.parmtop << _END_ > ptraj_CA_out2pdb1.out
ptraj prot_0wat.parmtop << _END_ > ptraj_CA_out2pdb1.out
trajin  prod_vac.nc 1 1
strip !@CA
trajout prod_CA.pdb pdb
go
_END_

