import sys
from prokka_helpers import *
from basic_helpers import *

prokka = "/home/fsepulveda/Documents/PROGRAMS/prokka-1.12/bin/prokka"

def simple_sequence(file):
	""" Toma un archivo fasta y devuelve la secuencia """
	with open(file) as f:
		lines = f.read()
		sequences = [line for line in lines.split('\n') if not ">" in line]
		sequence = ''.join(map(str,sequences))
	return sequence


cassette = simple_sequence(sys.argv[1])

print len(cassette)

execute_prokka(prokka, "OUT_TEST", sys.argv[1])

ffn, gff, fna, faa = prokka_files("OUT_TEST")

faa_dict = fasta2dict(faa)

ordered_keys = []
for key in faa_dict.keys():
	print key
	ordered_keys.append(key)

print

ordered_keys = sorted(ordered_keys)
for k in ordered_keys:
	sequence = get_sequence(faa_dict, k)
	print sequence