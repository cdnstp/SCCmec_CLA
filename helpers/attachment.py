import os
import re
import subprocess
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #

def attr_blast_search(blast_exe, database, query):
	""" Ejecuta BLAST usando database+string """
	outfmt = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"

	cmd = [blast_exe, "-outfmt", outfmt, "-word_size", "6", "-db", database]

	proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
								stdout=subprocess.PIPE,
								stderr=subprocess.STDOUT)
	
	results, err = proc.communicate(query.encode())

	return results, err

def core_finder(dna_sequence):
	""" returns dict with keys as attr-sequences and values as positional data """
	pass
	core = 'TATCATA[AT]'
	matches = re.finditer('{}'.format(core), dna_sequence)

	attr_location = {}
	for match in matches:
		span = match.span()
		start, end = span[0]-29, span[1]+23
		if (start >= 1000) and (start <= 100000):
			hypothetical_attr = dna_sequence[start:end]
			attr_location[hypothetical_attr] = [start, end]

	return attr_location

def attr_filter(attr_location, blastn_exe, attr_database_path):
	blast_hits = []
	for k, v in attr_location.items():
		out, err = attr_blast_search(blastn_exe, attr_database_path, k)
		if out:
			out = out.decode('utf8')
			hits = [line.split('\t') for line in out.split('\n')]
			hits = [line for line in hits if line != [""]]
			for hit in hits:
				qcov = (((float(hit[5])-float(hit[4]))+1)/float(hit[1]))*100
				if qcov > 30.0:
					hit.append(str(format(qcov, '.2f')))
					hit.append(v)
					blast_hits.append(hit)

	if blast_hits != []:
		return blast_hits
	else:
		return None

def select_attr(blast_hits):
	sorted_matches = sorted(blast_hits, key=lambda x: [float(x[12]), float(x[10])], reverse=True)
	print(sorted_matches[0])
	best_hit = sorted_matches[0]
	return best_hit[13][0], best_hit[13][1], sorted_matches[0]

def get_attchment(nombre, output, template_dna, attr_database_path, blastn_exe):
	os.chdir(output)
	attr_location = core_finder(template_dna)
	blast_hits = attr_filter(attr_location, blastn_exe, attr_database_path)
	if blast_hits is not None:
		start, end, hit = select_attr(blast_hits)
		attr_sequence = template_dna[start:end]
		write_fasta(nombre, start, end, attr_sequence)
		return start, end, hit

def write_fasta(nombre, start, end, attr_sequence):
    filename = 'attr_{}.fasta'.format(nombre)
    header = '>attr_{}_{}-{}'.format(nombre, str(start), str(end))
    with open(filename, "w") as f:
        f.write("{}\n{}\n".format(header, attr_sequence))
            # with open(filename, 'w') as f:
    #     f.write('{}\n{}\n'.format(header, attr_sequence))

