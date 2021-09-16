import sys

if len(sys.argv) != 2:
	print("Usage: ./simpleCount <path_to_bed_file>")

path = sys.argv[1]

tot = 0
with open(path, 'r') as f:
	for line in f:
		line = line.strip()
		line = line.split('\t')
		s = int(line[1])
		e = int(line[2])
		tot += e - s

print(tot)