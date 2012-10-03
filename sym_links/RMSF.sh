#!/bin/bash
# To find the rmsf per residue of an mdcrd or nc file
# Can try atomicfluct out RMSF.dat :1-486@CA byres
# http://ringo.ams.sunysb.edu/index.php/2012_AMBER_Tutorial_with_Biotin_and_Streptavidin#RMSD_Plots

ptraj ../prot_0wat.parmtop << end
trajin ../prod_vac.nc
atomic out RMSF.dat byres
end

