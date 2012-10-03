Created by Jonathan Teo 20/07/2012

Copy the whole folder into the simulation directory containing the pdb files
Change the Initial parameters in ala_scan.sh,
	including which residues to mutate,
	name of pdb file
	number of chains for splitpdb.py
Make sure this files are present in the parent folder
	The .nc trajectories file
	The master pdb file
	The pdb file with all water removed and tleaped (prot_0wat.pdb)

===> Run ala_scan.sh

Copy the ala_scan folder to a cluster

===> Run run.sh

