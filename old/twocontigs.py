# TEST SLIDING WINDOW 
from itertools import islice
import sys
import os
import subprocess

def fasta2dict(file):
	fastadict = {}
	with open(file) as file_one:
		for line in file_one:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				active_sequence_name = str(line.split()[0][1:])
				if active_sequence_name not in fastadict:
					fastadict[active_sequence_name] = []
				continue
			sequence = line
			fastadict[active_sequence_name].append(sequence)

	return fastadict

def get_sequence(sequences_dict, sequence_id):
	sequence = ""
	for i in sequences_dict.get(sequence_id)[:]:
		sequence += ''.join(i)

	return sequence

def window(seq, n):
	"Returns a sliding window (of width n) over data from the iterable"
	"   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
	it = iter(seq)
	result = tuple(islice(it, n))
	if len(result) == n:
		yield result
	for elem in it:
		result = result[1:] + (elem,)
		yield result


def blast(blast_exe, database, query):
	outfmt = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"
	cmd = [blast_exe, "-outfmt", outfmt, "-db", database, "-word_size", "8"]
	proc = subprocess.Popen(cmd,
		stdin=subprocess.PIPE,
		stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT)
	results, err = proc.communicate(query)

	return results, err

def makeblastdb(makeblastdb_exe, working_dir, input_seq, nucl, db_name):
	os.chdir(working_dir)
	proc = makeblastdb_exe+' -in '+input_seq+' -dbtype '+nucl+' -out '+db_name
	os.system(proc)

	return db_name

def create_dir(base, dir_name):
	new_dir = os.path.join(base, dir_name)
	if not os.path.isdir(new_dir):
		os.makedirs(new_dir)
	return new_dir


blast_exe = "/usr/bin/blastn"
makeblastdb_exe = "/usr/bin/makeblastdb"

raw_data = os.getcwd()

attr_database_dir = create_dir(raw_data, "attr_database")

attr_path = "/home/fsepulveda/Desktop/working_dir/input/attR_database.txt"
database = makeblastdb(makeblastdb_exe, attr_database_dir, attr_path, "nucl", "attr_nucl_db")

fna = sys.argv[1]
contigs_dict = fasta2dict(fna)

contig_id = "gnl|10|saureus_1"

template_dna = get_sequence(contigs_dict, contig_id)

attr = "CCGCATCGTTAAATGATACGCAGAGGCGTATCATAAGTAATGAGGTTCATGATTTTTG"
size_window = len(attr)
print("ZISE: ", size_window)


windows = ["".join(x) for x in window(template_dna, size_window)]
core = "TATCATAA"


for i in windows:
	if core in i:
		result, err = blast(blast_exe, database, i)
		print result
		for hit in result.split('\n'):
			if hit:
				args = [arg for arg in hit.split('\t')]
				porcentaje = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100
				if porcentaje >= 60.0:
					print hit


