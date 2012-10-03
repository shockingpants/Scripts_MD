#!/bin/bash
# ala.sh
# Prep files for alanine scan
workdir=`pwd`

###############################
##### User Input
###############################

if [ -z $1 ]; then
		echo "Jon: Please choose an option/argument!(-s -n)"
		exit 1
fi

while getopts ":ns" opt; do
  case $opt in
    n)
	  # Creates a new ala_scan folder
	  # Deletes any existing ala_scan
      rm -r ala_scan
	  cp -r /Volumes/HDD/teojy/Scripts/ala_scan .
      ;;
	s)
      # Copies over ala_scan scripts only
	  cd ala_scan
	  cp -r /Volumes/HDD/teojy/Scripts/ala_scan/ .
	  ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

