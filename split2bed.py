# Written by Daniel R. Olson 
#		  on Sept 7 2021


#Usage: split2bed <input_dir> <-trf|-ultra> <input_file_extension> <output_name> 

# This will split a fasta input_file into several files, where each file is
# a subsequence of length <symbols_per_file>


import sys
import os
from os import listdir
from os.path import isfile, join

class SplitFileReader():

	def __init__(self, ultra, seq_leng, seq_overlap):
		self.ultra = ultra
		self.seqLeng = seq_leng
		self.seqOverlap = seq_overlap

		self.curSeq = ''
		self.curPos = 0
		self.seqs = dict()
		self.seqNames = []

	def outputToFile(self, filePath):
		with open(filePath, 'w') as file:
			for seq in self.seqNames:
				for line in self.seqs[seq]:
					print(seq, line[0], line[1], line[2], sep='\t', file=file)

	def readULTRAFile(self, filePath):
		seqName = ''
		cpos = self.curPos - self.seqOverlap
		if cpos < 0:
			cpos = 0
		stage = 0

		with open(filePath, 'r') as file:
			start = 0
			end = 0
			score = 0
			for line in file:
				line = line.strip()
				line = line.split(' ')

				if stage == 0:
					if line[0] == '"SequenceName":':
						stage = 1
						seqName = line[1][1:-2]
						if seqName != self.curSeq:
							self.curPos = 0
							cpos = 0
							self.curSeq = seqName
							self.seqNames.append(seqName)
							self.seqs[seqName] = []
			    #looking for start
				elif stage == 1:
					#print('stage1', line[0])
					if line[0] == '"Start":':
						stage = 2
						start = int(line[1][:-1])
						#print('stage1b')
			    #looking for end
				elif stage == 2:
					#print('stage 2')
					if line[0] == '"Length":':
						stage = 3
						end = start + int(line[1][:-1])
						#print('stage2b')
				#looking for score
				elif stage == 3:
					#print('stage 3')
					if line[0] == '"Score":':
						stage = 1
						score = float(line[1][:-1])
						self.seqs[seqName].append((start + cpos, end + cpos, score))
						#print("stage3b")




	def readTRFFile(self, filePath):
		seqName = ''
		cpos = self.curPos - self.seqOverlap
		if cpos < 0:
			cpos = 0
		stage = 0

		with open(filePath, 'r') as file:

			for line in file:
				line = line.strip()
				line = line.split(' ')
				if stage == 0:
					if line[0] == 'Sequence:':
						stage = 1
						seqName = line[1]
	
						if seqName != self.curSeq:
							self.curPos = 0
							cpos = 0
							self.curSeq = seqName
							self.seqNames.append(seqName)
							self.seqs[seqName] = []
	
				elif len(line) > 10:
					lineStart = int(line[0]) - 1
					lineEnd = int(line[1]) - 1
					lineScore = float(line[7])
					self.seqs[seqName].append((lineStart + cpos, lineEnd + cpos, lineScore))

				else:
					if len(line) > 0:
						if line[0] == 'Sequence:':
							stage = 1
							seqName = line[1]
		
							if seqName != self.curSeq:
								self.curPos = 0
								cpos = 0
								self.curSeq = seqName
								self.seqNames.append(seqName)
								self.seqs[seqName] = []


	def moveForMissingFile(self):
		self.curPos += self.seqLeng

	def readFile(self, filePath):
		if self.ultra:
			self.readULTRAFile(filePath)
		else:
			self.readTRFFile(filePath)
		self.curPos += self.seqLeng


if len(sys.argv) != 7:
	print("Usage: split2bed <input_dir> <-trf|-ultra> <input_file_extension> <output_name> <seq_leng> <seq_overlap>")
	print("Example: split2bed ./split_files/ -trf dat trf_beds 50000 6000")
	exit()

input_dir = sys.argv[1]
input_data = sys.argv[2]
extension = sys.argv[3]
output_name = sys.argv[4]
seq_leng = int(sys.argv[5])
seq_overlap = int(sys.argv[6])

if input_data != '-trf' and input_data != '-ultra':
	print("Usage: split2bed <input_dir> <-trf|-ultra> <input_file_extension> <output_name> <seq_leng> <seq_overlap>")
	print("Example: split2bed ./split_files/ -trf dat trf_beds 50000 6000")
	exit()


files = [file for file in listdir(input_dir) if isfile(join(input_dir, file)) and file.endswith(extension)]

numToFile = dict()
max_file = 0

for f in files:
	file = join(input_dir, f)
	num = int(f[f.find('_') + 1:f.find('.')])
	if num > max_file:
		max_file = num
	numToFile[num] = file

ultra = True
if input_data == '-trf':
	ultra=False

reader = SplitFileReader(ultra, seq_leng, seq_overlap)
for i in range(max_file):
	if i not in numToFile:
		reader.moveForMissingFile()
	else:
		reader.readFile(numToFile[i])
		
reader.outputToFile(output_name)
