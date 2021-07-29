# Written by Daniel R. Olson 
#		  on Thu Jul 29 2021


#Usage: readtimes <dir_path> <file extension>

# reads the output of '/usr/bin/time -v'  output in a directory and
# then outputs the results in csv format to stdev


import sys
import os
from os import listdir
from os.path import isfile, join

if len(sys.argv) != 3:
	print("Usage: readtimes <dir_path> <file extension>")
	print("Example: readtimes ~/documents/ txt")
	exit()



