#!/bin/bash
# extract_avg_dG.sh
# Extracts the average binding energy (and possibly other average energies) for alascan and do a comparison

rm -r ala_energy 

##################################
### File with original energy data
##################################
mkdir ala_energy
cd ala_energy
ener_dir=$(pwd)
cd ..
ala_dir=$(pwd)
# Checks and creates list of directories
if [[ ! -f folder_list.txt ]]; then
	ls -1 | grep -P '\d+_\w' >> folder_list.txt
	cat folder_list.txt | sort -n > fold.txt
	mv fold.txt folder_list.txt
else
	cat folder_list.txt | sort -n > fold.txt
	mv fold.txt folder_list.txt
fi
list_dir=$(echo $ala_dir/folder_list.txt)

echo Tot_energy STD STDE_mean > "$ener_dir/wt_totalbind.dat"
grep 'DELTA G binding' ../mmgbsa/03_MMPBSA_binding/FINAL_RESULTS_MMPBSA.dat | awk '{print $5, $7, $8}' >> "$ener_dir/wt_totalbind.dat"


###############################################
### Extract energy values from the ala scan files
###############################################

echo resn Tot_energy STD STDE_mean > "$ener_dir/mut_totalbind.dat"
echo resn ddG > "$ener_dir/mut_ddG.dat"


for i in $(cat folder_list.txt); do
	cd ${i}
	if [[ ! -f FINAL_RESULTS_MMPBSA.dat ]]; then
		echo $i >> $ener_dir/missing.log
		cd ..
		continue
	fi
	res=$(echo ${i} | sed 's|/||g')
# To account for difference between amber11 and amber12 (mpi)
# The one with mutant is from amber 11
	grep 'DELTA DELTA G binding'  FINAL_RESULTS_MMPBSA.dat | awk '{ if ($6 == "MUTANT:)") print "'$res'", $12*-1; else print "'$res'", $11 }' >> "$ener_dir/mut_ddG.dat"
	awk 'NR > 95 {print $0}' FINAL_RESULTS_MMPBSA.dat |grep '^ DELTA G binding'| awk '{print "'$res'", $5,$7,$8}' >> "$ener_dir/mut_totalbind.dat"
	cd ..
done

# Get number of data

hi=$(cat "$ener_dir/mut_ddG.dat" | wc -l)
hi=$(($hi-1))
echo $hi > $ener_dir/number_of_success.dat

##############################################
### Plot the energy in terms of a bar graph
#############################################
#par=$(pwd|sed 's|.*P53/\(.*\)/ala_scan|\1|')
cd ../ ; par=$(basename `pwd`) ; cd - > /dev/null
cd ala_energy
python2.7 ../plot_ala.py $par





