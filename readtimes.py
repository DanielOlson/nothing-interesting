# Written by Daniel R. Olson 
#		  on Thu Jul 29 2021


#Usage: readtimes <dir_path> <file extension>

# reads the output of '/usr/bin/time -v'  output in a directory and
# then outputs the results in csv format to stdev
# This is not a fancy program,
#	don't use it on anything but very sanitized files


import sys
import os
from os import listdir
from os.path import isfile, join

if len(sys.argv) != 3:
	print("Usage: readtimes <dir_path> <file extension>")
	print("Example: readtimes ~/documents/ .txt")
	exit()

in_dir = sys.argv[1]
extension = sys.argv[2]

all_files = [join(in_dir, file) for file in listdir(in_dir) if isfile(join(in_dir, file)) and file.endswith(extension)]


class Time():
	def __init__(self, fileName="<undefined>"):
		self.fileName = fileName
		self.user_time = ""
		self.system_time = ""
		self.pct_cpu = ""
		self.elapsed_time = ""
	#	self.avg_sharedtextsize = ""
		self.max_ressetsize = ""
		self.exit_status = ""

	def headerString():
		return "FileName\tUserTime\tSystemTime\tPctCPU\tWallTime\tMaxResSetSize\tExitStatus"

	def __str__(self):
		return "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(self.fileName,
			self.user_time, 
			self.system_time, 
			self.pct_cpu, 
			self.elapsed_time, 
			self.max_ressetsize,
			self.exit_status)

	def readLine(self, line):
		if "\n" in line:
			line = line[:-1]
		if "User time" in line:
			val = line[line.find(': ') + 2:]
			self.user_time = val

		if "System time" in line:
			val = line[line.find(': ') + 2:]
			self.system_time = val

		if "Percent of CPU" in line:
			val = line[line.find(': ') + 2:]
			self.pct_cpu = val	

		if "Elapsed (wall clock) time" in line:
			val = line[line.find(': ') + 2:]
			self.elapsed_time = val

		if "Maximum resident set size" in line:
			val = line[line.find(': ') + 2:]
			self.max_ressetsize = val

		if "Exit status" in line:
			val = line[line.find(': ') + 2:]
			self.exit_status = val

times = []
for f in all_files:
	newTime = Time(f)
	with open(f, 'r') as file:

		for line in file:
			newTime.readLine(line)
		times.append(newTime)


print(Time.headerString())
for t in times:
	print(str(t))



