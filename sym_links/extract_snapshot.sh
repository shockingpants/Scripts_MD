#!/bin/bash
# Extracts a snapshot from a netcdf (.nc) file


echo "To check total number of frames..."
ptraj prot_0wat.parmtop << end
trajin prod_vac.nc 1 1
end


read -p "Which frame would you like to extract?  " frame
ptraj prot_0wat.parmtop << end
trajin prod_vac.nc $frame $frame
trajout "prot_$frame.pdb" pdb
end
