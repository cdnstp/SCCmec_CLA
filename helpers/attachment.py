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
	
	results, err = proc.communicate(query)

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
	for k, v in attr_location.iteritems():
		out, err = attr_blast_search(blastn_exe, attr_database_path, k)
		if out:
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

def write_fasta(nombre, start, end, attr_sequence):
	filename = 'attr_{}.fasta'.format(nombre)
	header = '>attr_{}_{}-{}'.format(nombre, str(start), str(end))
 	with open(filename, 'w') as f:
 		f.write('{}\n'.format(header))
 		f.write('{}\n'.format(attr_sequence))

def get_attchment(nombre, output, template_dna, attr_database_path, blastn_exe):
	os.chdir(output)
	attr_location = core_finder(template_dna)
	blast_hits = attr_filter(attr_location, blastn_exe, attr_database_path)
	if blast_hits is not None:
		start, end, hit = select_attr(blast_hits)
		attr_sequence = template_dna[start:end]
		write_fasta(nombre, start, end, attr_sequence)
		return start, end, hit


# def get_attchment(nombre, output, template_dna, attr_database_path, blastn_exe):
# 	os.chdir(output)
# 	core = "TATCATAA"
# 	matches = re.finditer("{}".format(core), template_dna)
# 	print(matches)
# 	print(type(matches))
# 	salida = ""
# 	atts_location = {}
# 	for match in matches:
# 		span = match.span()
# 		print(span)
# 		desde = span[0]-28
# 		hasta = span[1]+22

# 		# Para evitar considerar attL2, posicion de attr encontrado > 1000pb
# 		if hasta >= 1000 and hasta <= 100000:

# 			hypothetical_attR = template_dna[desde:hasta]

# 			atts_location[hypothetical_attR] = [desde, hasta]

# 			out, err = attR_BLAST(blastn_exe, attr_database_path, hypothetical_attR)

# 			row = [s.split('\t') for s in out.split('\n')]

# 			lines = [x for x in row if x != [""]]

# 			for line in lines:

# 				args = [arg for arg in line]

# 				coverage = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100

# 				if coverage >= 30.0:
# 					args.append(str(format(coverage, '.2f')))
# 					args.append(str(desde))
# 					args.append(str(hasta))
# 					args[0] = hypothetical_attR
# 					salida += "\t".join([x for x in args])+'\n'

# # --------------------------------------------------------------------------- #
# #   En caso de encontrar solo 1 posible attr
# 	if len(atts_location) == 1:
# 		print("1 attR encontrado")
# 		attr_sequence = atts_location.keys()[0]
# 		attR_inicio = atts_location.values()[0][0]
# 		attR_final = atts_location.values()[0][1]
# 		hit = "attR_"+nombre+"_"+str(attR_inicio)+"-"+str(attR_final)
# 		with open("attR_"+nombre+".fasta", "w") as f:
# 			f.write(">"+hit+"\n")
# 			f.write(attr_sequence+"\n")

# 		return int(attR_inicio), int(attR_final), hit
# # --------------------------------------------------------------------------- #
# # --------------------------------------------------------------------------- #
# #   En caso de haber mas de 1 posible attr, filtrar por mayor coverage
# 	if salida:
# 		lines = [s.split('\t') for s in salida.split('\n')]
# 		list2 = [x for x in lines if x != [""]]
# 		print(list2)
# 		list2.sort(key=lambda x: float(x[12]), reverse=True)
# 		print("Best Hit attR: ", list2[0])
# 		attr_sequence = list2[0][0]

# 		print attr_sequence
# 		attR_inicio = list2[0][13]
# 		attR_final = list2[0][14]
# 		hit = "attR_"+nombre+"_"+str(attR_inicio)+"-"+str(attR_final)
# 		with open("attR_"+nombre+".fasta", "w") as f:
# 			f.write(">"+hit+"\n")
# 			f.write(attr_sequence+"\n")
# 		return int(attR_inicio), int(attR_final), hit

# 	else:
# 		print("attR not found")
# 		return None
