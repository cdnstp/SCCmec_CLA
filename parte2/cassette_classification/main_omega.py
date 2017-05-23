import sys
import os
import itertools
import subprocess

path_to_mash = "/home/fsepulveda/Documents/parte2/cassette_classification/ex_programs/mash-Linux64-v1.1.1/mash"

def mash(path_to_mash, cassette1, cassette2, kmer):
	cmd = path_to_mash+" dist "+cassette1+" "+cassette2+" -k {0}".format(kmer)
	p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	output = p.stdout.read()
	print output
	return output



def main():
	path_to_files = os.listdir(sys.argv[1])

	output = []
	valor = "13"
	salida_base = "/home/fsepulveda/Documents/parte2/cassette_classification/"

	print("\n"+"\n"+"\t"+"\t"+"\t"+valor+"\n"+"\n"+"-"*78+"\n"+"-"*78)

	for par in itertools.product(path_to_files, repeat=2):

		p1 = par[0]
		p2 = par[1]
		
		out_mash = mash(path_to_mash, p1, p2, valor)
		output.append(out_mash+'\n')

	outfile = "raw_mash_output_k{0}_s1000.txt".format(valor)
	out_path = "".join([salida_base, outfile])
	with open(out_path, "w") as f:
		for i in output:
			f.write(i)

###
###
"""
./cassette_classification/ex_programs/mash-Linux64-v1.1.1/mash sketch -l ./cassette_classification/path_to_cassettes.txt -o chunk_s5000_k19 -s 5000 -k 19
"""
###
###

if __name__ == '__main__':
	main()
