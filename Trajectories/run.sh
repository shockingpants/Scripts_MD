#!/bin/bash
# run.sh
# Change variables to the relevant directories and run file

# Change _DIR_, _PWD_ in all files that need changing

santa_nam=${PWD##*/}
sed -i "s|_DIR_|`pwd`|" gpu.sh
sed -i "s|_DIR_|`pwd`|" comp_mdcrd.sh
sed -i "s|_NAME_|$santa_nam|" gpu.sh
sed -i "s|_NAME_|$santa_nam|" comp_mdcrd.sh


#cd mmg*/03*
#sed -i "s|_DIR_|`pwd`|" 03.sh
#sed -i "s|_NAME_|$santa_nam|" 03.sh
#
#cd ../..

chmod 777 gpu.sh
qsub gpu.sh
