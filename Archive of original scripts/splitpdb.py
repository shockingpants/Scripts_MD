#!/opt/local/bin/python2.7

'''
Yun Liu
6/14/2010

Splits pdb file into targets, based on TER
WARNING: overwrites the targets!!!

usage:
./splitpdb.py rec.pdb lig.pdb < re.pdb
./splitpdb.py 1.pdb 2.pdb 3.pdb 4.pdb [etc...] < big_protien.pdb
'''

import sys
# Jon added this
from strtup import *

# get targets
targets = sys.argv[1:] # sys.argv reads in the arguments placed into command line. The .py file is considered index 0 and counts as an input. so typically
		       # arguments start from index one. In this case, everything after .py. 

#===============================
#      Extracting file
#===============================
print 'Extracting file content...'
i=0
output = [[]]
for line in sys.stdin:
	# ignore remark line

	if not line.startswith("REMARK"):
		print i
		output[i].append(line)
		print output
		if line.startswith("TER"):
			output.append([])
			i += 1
		elif line.startswith("END"):
			break # end		
		
# clear empty chains
to_delete = []
for i in range(len(output)):
	if len(output[i]) == 0:
		to_delete.append(i)
for i in to_delete[::-1]: # delete in reverse
	del output[i]
	

len_t,len_o = len(targets), len(output)
if len_t == 0 and len_o != 0:
	sys.stderr.write("WARNING: no targets specified; writing to TARGETX.pdb\n")
	for i in range(len_o):
		targets.append("TARGET%d.pdb" % (i+1))
	len_t = len_o

if len_t != len_o:
	sys.stderr.write("ERROR: %d targets does not match %d chains in pdb\n" % (len_t,len_o))
	sys.exit()

for i in range(len_t):
	f = open(targets[i],'w')
	for line in output[i]:
		f.write(line)
	f.close()

#sys.stderr.write("LOG: Wrote %d targets, files written were: %s\n" %(len_t,str(targets)))
