#!/bin/bash
# compress.sh
# Compress all mccrd files
if [[ -z $1 ]]; then
	echo 'Please choose a compression rate. -# [1-9(slowest)]'
	exit
fi

echo 'Compressing... '
for j in $(ls -1|grep '.*\.mdcrd$'); do
	# Checks for existence of .gz
	if [ -f ${j} -a ! -f ${j}.gz ]; then
		echo $j `date`>> compress.log
		echo "Compressing ${j}"
		gzip $1 ${j} &
	fi
done
