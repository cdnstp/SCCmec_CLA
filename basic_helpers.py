# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #

def simple_sequence(file):
	""" Toma un archivo fasta y devuelve la secuencia """
	with open(file) as f:
		lines = f.read()
		sequences = [line for line in lines.split('\n') if not ">" in line]
		sequence = ''.join(map(str,sequences))
	return sequence


# ------------------------------------------------------------------------- #

def reverse_complement(sequence):
	""" Transforma la secuencia nucl en su reverso complementario """
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	bases = list(sequence)
	bases = [complement[base] for base in bases]
	complement = ''.join(bases)
	reverse = complement[::-1]
	return reverse


# ------------------------------------------------------------------------- #

def create_dir(base_path, dir_name):
	""" Crea un directorio y devuelve fullpath """
	import os
	new_dir = os.path.join(base_path, dir_name)
	if not os.path.isdir(new_dir):
		os.makedirs(new_dir)
	return new_dir


# ------------------------------------------------------------------------- #

def fasta2dict(file):
	""" Transforma archivos multi-fasta en un dict de python """
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


# ------------------------------------------------------------------------- #

def get_sequence(sequences_dict, sequence_id):
	""" Obtiene secuencia desde dict-seq usando la key """
	sequence = ""
	for i in sequences_dict.get(sequence_id)[:]:
		sequence += ''.join(i)
	return sequence


# ------------------------------------------------------------------------- #

def ordenar_contigs(orfx, mec, ccr):
	""" Si core-elements en diferentes contigs, saber cual es el primero (left, orfx)"""
	ids = [orfx, mec, ccr]
	print len(set(ids))
	if not len(set(ids)) > 2:
		if orfx == mec:
			contig_one = orfx
			contig_two = ccr
			print("orfx + mec")
			return contig_one, contig_two
		if orfx == ccr:
			contig_one = orfx
			contig_two = mec
			print("orfx + ccr")
			return contig_one, contig_two
		else:
			contig_one = orfx
			contig_two = mec
			print("mec + ccr")
			return contig_one, contig_two
	else:
		print("Cannot be sorted")
		sys.exit()


# ------------------------------------------------------------------------- #

def checkContig(lst):
	""" Comparar lista (de contigs) para saber si son los mismos """
	return lst[1:] == lst[:-1]


# ------------------------------------------------------------------------- #



def checkSense(sequence, contig):
	import re
	""" Return always attL-attR left to right """
	if re.findall("{pattern}".format(pattern=sequence), contig):
		return contig

	contig = reverse_complement(contig)
	if re.findall("{pattern}".format(pattern=sequence), contig):
		return contig


# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #

def get_orfX_pos(sequence, query, name):
	import re
	if re.findall("{pattern}".format(pattern=query), sequence):
		for i in re.findall("{pattern}".format(pattern=query), sequence):
			x = "(+) %s located_at:" % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]
			texto = "(+) {0} located_at:".format(name)
			inicio = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)][0][0]
			final = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)][0][1]
			#print texto, inicio, final
			return texto, inicio, final
	else:
		sys.exit("Error")


# ------------------------------------------------------------------------- #

def sequence_position(sequence, query, name):
	import re
	sentido = re.search(query, sequence)
	if sentido is not None:
		texto = "(+) {0} located_at:".format(name)
		inicio, final = sentido.span()[0], sentido.span()[1] 
		print texto, inicio, final
		return texto, inicio, final

	antisentido = re.search(reverse_complement(query), sequence)
	if antisentido is not None:
		texto = "(-) {0} located_at:".format(name)
		inicio, final = antisentido.span()[1], antisentido.span()[0]
		print texto, inicio, final
		return texto, inicio, final
	else:
		print(name, "Not Found")
		import sys
		sys.exit()


# ------------------------------------------------------------------------- #