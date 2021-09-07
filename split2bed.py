# Written by Daniel R. Olson 
#		  on Sept 7 2021


#Usage: split2bed <input_dir> <-trf|-ultra> <input_file_extension> <output_name> 

# This will split a fasta input_file into several files, where each file is
# a subsequence of length <symbols_per_file>


import sys
import os
from os import listdir
from os.path import isfile, join

if len(sys.argv != 5):
	print("Usage: split2bed <input_dir> <-trf|-ultra> <input_file_extension> <output_name> ")
	print("Example: split2bed ./split_files/ -trf dat trf_beds")
	exit()

input_dir = sys.argv[1]
input_data = sys.argv[2]
extension = sys.argv[3]
output_name = sys.argv[4]

if input_data != '-trf' and input_data != '-ultra':
	print("Usage: split2bed <input_dir> <-trf|-ultra> <input_file_extension> <output_name> ")
	print("Example: split2bed ./split_files/ -trf dat trf_beds")
	exit()


files = [join(input_dir, file) for file in listdir(input_dir) if isfile(join(input_dir, file)) and file.endswith(extension)]

for f in files:
	print(f)
