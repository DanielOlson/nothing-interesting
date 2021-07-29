# Written by Daniel R. Olson 
#		  on July 22 2021


#Usage: splitit <input_file> <output_name> <symbols_per_file>

# This will split a fasta input_file into several files, where each file is
# a subsequence of length <symbols_per_file>


import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from math import ceil



if len(sys.argv) != 5:
	print("Usage: splitit <input_file> <output_name> <symbols_per_file> <overlap>")
	print("Example: splitit input.fa genome 1000")
	exit()


input_file = sys.argv[1]
output_name = sys.argv[2]
seqLength = int(sys.argv[3])
overlap = int(sys.argv[4])

#Read the input file
records = list(SeqIO.parse(input_file, "fasta"))

files_counter = 0

for r in records:
	name = r.id
	seq = r.seq
	bins = ceil(len(seq) / seqLength)

	for i in range(bins):
		stop = False
		start = i * seqLength
		if i > 0:
			start = start - overlap
		end = (i + 1) * seqLength
		end += overlap
		if end > len(seq):
			end = len(seq)
			stop = True

		newrec = SeqRecord(Seq(seq[start:end]), id=r.id, name=r.name, description=r.description)

		with open(output_name + '_' + str(files_counter) + '.fa', 'w') as output_handle:
			SeqIO.write(newrec, output_handle, 'fasta')
		files_counter += 1

		if stop == True:
			break

	print(bins)