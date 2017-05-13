import os
import sys
import regex as re
import subprocess

from itertools import islice

def config():
	import ConfigParser
	config = ConfigParser.ConfigParser()
	config.readfp(open(r'config.txt'))
	prokka_exe = config.get('configuration', 'prokka')
	blastn_exe = config.get('configuration', 'blastn')
	blastp_exe = config.get('configuration', 'blastp')
	makeblastdb_exe = config.get('configuration', 'makeblastdb')
	inputFiles = config.get('configuration', 'input')
	seqFile = config.get('configuration', 'seqFile')
	ccr = config.get('configuration', 'ccr')
	orfX = config.get('configuration', 'orfX')
	mecA = config.get('configuration', 'mecA')
	attr_db = config.get('configuration', 'attr_db')

	return prokka_exe, blastn_exe, blastp_exe, makeblastdb_exe, inputFiles, seqFile, ccr, orfX, mecA, attr_db


def simple_sequence(file):
    with open(file) as f:
        lines = f.read()
        sequences = [line for line in lines.split('\n') if not ">" in line]
    sequence = ''.join(map(str,sequences))

    return sequence

def blastAlign(blast_exe, query, subject):
	formato = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"
	process = subprocess.Popen([blast_exe, "-word_size", "7", "-query", query, "-subject", subject, "-outfmt", formato], stdin=subprocess.PIPE,
		stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT)
	out, err = process.communicate()
	#print out
	salida = ""
	for hit in out.split('\n'):
		if hit:
			args = [arg for arg in hit.split('\t')]
			#print args
			porcentaje = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100
			#print porcentaje
			""" OJOOOOOO  """
			if porcentaje >= 60.0:
				salida += hit

	return salida, err

def blastAlign2(blast_exe, query, subject):
	formato = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"
	process = subprocess.Popen([blast_exe, "-word_size", "6", "-query", query, "-subject", subject, "-outfmt", formato], stdin=subprocess.PIPE,
		stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT)
	out, err = process.communicate()

	salida = ""
	for hit in out.split('\n'):
		if hit:
			args = [arg for arg in hit.split('\t')]
			#print args
			porcentaje = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100
			#print porcentaje
			""" OJOOOOOO  """
			if porcentaje >= 20.0:
				salida += hit+'\t'+str(format(porcentaje, '.2f'))+"\n" 

	return salida, err

def blast(blast_exe, database, query):
	outfmt = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"
	cmd = [blast_exe, "-outfmt", outfmt, "-db", database, "-word_size", "7"]
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


def prokka_files(location):
	for file in os.listdir(location):
		if file.endswith(".ffn"):
			ffn = os.path.join(location, file)
		if file.endswith(".gff"):
			gff = os.path.join(location, file)
		if file.endswith(".fna"):
			fna = os.path.join(location, file)
		if file.endswith(".faa"):
			faa = os.path.join(location, file)

	return ffn, gff, fna, faa


def execute_prokka(prokka_exe, output_prokka, contigs):
	kingdom = "Bacteria"
	genus = "Staphylococcus"
	locustag = "saureus"
	cmd_prokka = prokka_exe+" --kingdom "+kingdom+ \
					" --outdir "+output_prokka+ \
					" --quiet "+ \
					" --genus "+genus+ \
					" --locustag "+locustag+ \
					" --centre 10 --compliant"+ \
					' '+contigs
	os.system(cmd_prokka)


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


#import string
#import random
#def id_generator(size=6, chars=string.ascii_lowercase + string.digits):
#	return ''.join(random.choice(chars) for _ in range(size))


def reverse_complement(seq):
	"""takes a sequence to get the reverse complement of it"""
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	bases = list(seq)
	bases = [complement[base] for base in bases]
	complement = ''.join(bases)
	reverse = complement[::-1]
	return reverse

def get_contig(gff, query):
	with open(gff) as f:
		header = []
		lines = f.readlines()
		for line in lines:
			if query in line:
				header.append(line.split()[0])

	return header[0]


def checkContig(lst):
	return lst[1:] == lst[:-1]


def checkSense(sequence, contig):
	"""Return always attL-attR left to right"""
	if re.findall("{pattern}".format(pattern=sequence), contig):
		return contig
	contig = reverse_complement(contig)
	if re.findall("{pattern}".format(pattern=sequence), contig):
		return contig


def create_dir(base, dir_name):
	new_dir = os.path.join(base, dir_name)
	if not os.path.isdir(new_dir):
		os.makedirs(new_dir)
	return new_dir


def find_position(sequence, query, name):
	if re.findall("{pattern}".format(pattern=query), sequence):
		for i in re.findall("{pattern}".format(pattern=query), sequence):
			x = "(+) %s located_at:" % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]
			texto = "(+) {0} located_at:".format(name)
			inicio = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)][0][0]
			final = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)][0][1]
			#print texto, inicio, final
			return texto, inicio, final
	else:
		for i in re.findall("{pattern}".format(pattern=reverse_complement(query)), sequence):
			x = "(-) %s located_at:" % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]
			texto = "(-) {0} located_at:".format(name)
			inicio = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)][0][1]
			final = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)][0][0]
			#print texto, inicio, final
			return texto, inicio, final

# ------------------------------------------------------------------------- #

def findAtts(nombre, output, template_dna, att_actual_orfx, attr_database_path, blastn_exe):
	os.chdir(output)
# ------------------------------------------------------------------------- #
#      ?att; patron basico corresponde a los nucl finales del gen orfX      #
	
	atts_location = {}

	for i in re.findall("({s}){{s<=5}}".format(s=att_actual_orfx), template_dna):
		atts_location[template_dna[m.start(0):(m.end(0))]] = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), template_dna)]

	core = "TATCATAA"
	adjacent_core = "GA[ACTG]G"
	""" eliminar att encontrados que no posean la secuencia central de att"""
	for key in atts_location.keys():
		if not core in key:
			del atts_location[key]
	for k in atts_location.keys():
		if not re.findall("{s}".format(s=adjacent_core), k):
			del atts_location[k]

	""" ordenar elementos en att_location segun posicion (creciente)"""
	sites = sorted(atts_location.items(), key=lambda e: e[1][0])
	print sites


# ------------------------------------------------------------------------- #
#              Separar attL attR segun posicion                             #

	hipotetical_attL = {}
	hipotetical_attR = {}

	value_attl = int([x+y for (x,y) in sites[0][1]][0])

	for s in sites:
		print s
		fr = s[1][0][0]
		to = s[1][0][1]
		if fr+to > (value_attl+300):
			hipotetical_attR[s[0]] = s[1][0]
		else:
			hipotetical_attL[s[0]] = s[1][0]

		attL_list = sorted(hipotetical_attL.items(), key=lambda e: e[1][0])
		attR_list = sorted(hipotetical_attR.items(), key=lambda e: e[1][0])
		attL_sequence = template_dna[attL_list[0][1][0]-20:attL_list[0][1][1]+20]

		with open("attL_"+nombre+".fasta", "w") as f:
			f.write(">attL_"+nombre+"\n")
			f.write(attL_sequence+'\n')

		output_attR = "attR_%s.fasta" % nombre
		with open(output_attR, "w") as f:
			for attr in attR_list:
				f.write(">attR_"+nombre+"_"+str(attr[1][0]-20)+"-"+str(attr[1][1]+20)+"\n")
				f.write(str(template_dna[attr[1][0]-20:attr[1][1]+20])+'\n')

# ------------------------------------------------------------------------- #
#     Filtrar attachment site sequence Right (attR) basado en similitud     #

	attR_list_path = os.path.join(output, output_attR)
	print attR_list_path
	if not os.stat(attR_list_path).st_size == 0:
		result, err = blastAlign(blastn_exe, attR_list_path, attr_database_path)
		hit = result.split()[0]
		if hit:
			print("HIT: ", hit)
			attR_start = hit.split("_")[-1].split("-")[0]
			attR_end = hit.split("_")[-1].split("-")[1]
			coordinates = [attL_list[0][1][0], int(attR_end)]
			print "attL starting at: ", attL_list[0][1][0]
			print "attR ends at: ", attR_end
			print coordinates

		return coordinates, hit
	else:
		print("attR no encontrado") 

# ------------------------------------------------------------------------- #

def blast_parser(blast_exe, db, query, qname):
	blastfile, err = blast(blast_exe, db, query)

	salida = ""
	for hit in blastfile.split('\n'):
		if hit:
			args = [arg for arg in hit.split('\t')]
			porcentaje = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100
			if porcentaje >= 70.0:
				#print args
				#if float(args[10]) >= 60.0:
				print("Porcentaje de Alineamiento: ", porcentaje)
				print("Porcentaje de Identidad: ", args[10])
				print hit
				salida += hit
	if not salida:
		print("Query {0} not found in SeqFile").format(qname)
		return None
	else:
		bestHit = salida.split()[2]
		return bestHit

def ordenar_contigs(orfx, mec, ccr):
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

def attR(contig_sequence, raw_data, att_actual_orfx, nombre, blastn_exe, attr_database_path, sense):

	os.chdir(raw_data)
	atts_location = {}
	for i in re.findall("({s}){{s<=5}}".format(s=att_actual_orfx), contig_sequence):
		atts_location[contig_sequence[m.start(0):(m.end(0))]] = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), contig_sequence)]

	core = "TATCATAA"
	adjacent_core = "GA[ACTG]G"
	for key in atts_location.keys():
		if not core in key:
			del atts_location[key]
	for k in atts_location.keys():
		if not re.findall("{s}".format(s=adjacent_core), k):
			del atts_location[k]

	sites = sorted(atts_location.items(), key=lambda e: e[1][0])

	print("HIPOTETICAL SITES: ", sites)

	if len(sites) > 0:
		output_attR = "attR_%s_%s.fasta" % (nombre, sense)
		with open(output_attR, "w") as f:
			for attr in sites:
				print("attr: ", attr)
				attR_start = attr[1][0][0]
				print("attR_start: ", attR_start)
				attR_end = attr[1][0][1]
				print("attR_end: ", attR_end)

				f.write(">attR_"+nombre+"_"+str(attR_start-20)+"-"+str(attR_end+20)+"\n")
				f.write(str(contig_sequence[attR_start-20:attR_end+20])+'\n')

		print("PROBABLI ATTR")		
		print(str(contig_sequence[attR_start-20:attR_end+20])+'\n')

		attR_list_path = os.path.join(raw_data, output_attR)

		#os.system("cat {0}".format(attR_list_path))

		result, err = blastAlign2(blastn_exe, attR_list_path, attr_database_path)
		print type(result)
		
		if result:
			print result
			lines = [s.split('\t') for s in result.split('\n')]
			list2 = [x for x in lines if x != [""]]
			#print len(list2)
			from operator import itemgetter
			list3 = sorted(list2, key=itemgetter(12), reverse = True)
			hit = result.split()[0]
			print list3[0]
			hit = list3[0][0]
			print("HIT: ", hit)
			attR_start = hit.split("_")[-1].split("-")[0]
			attR_end = hit.split("_")[-1].split("-")[1]
			print "attR starts at: ", attR_start
			print "attR ends at: ", attR_end
			sccmec_rightEnd = contig_sequence[0:int(attR_end)]
			print("Largo sccmec right end: ", len(sccmec_rightEnd))
			#print hit
			return sccmec_rightEnd, hit 
		else:
			return "empty"
	else:
		return "empty"

# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
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
def get_score(word1, word2):
    shared_chars = [char for char in word1.lower() if char in word2.lower()]
    return len(shared_chars)




def main():
# ------------------------------------------------------------------------- #
#                      CONFIGURATION, OUTPUT SET UP                         #

	prokka_exe, blastn_exe, blastp_exe, makeblastdb_exe, inputFiles, seqFile, ccr, orfX, mecA, attr_db = config()

	orfx_base = simple_sequence(os.path.join(inputFiles, orfX))
	mecA_base = simple_sequence(os.path.join(inputFiles, mecA))
	ccr_base = simple_sequence(os.path.join(inputFiles, ccr))

	attr_database_path = os.path.join(inputFiles, attr_db)

	nombre = seqFile.split(".")[0]

	contigs = os.path.join(inputFiles, seqFile)

	working_dir = os.getcwd()

	output = "output_"+nombre

	raw_data = os.path.join(working_dir, output)

	output_prokka = os.path.join(raw_data, "output_prokka")

	execute_prokka(prokka_exe, output_prokka, contigs)

	ffn, gff, fna, faa = prokka_files(output_prokka)

# ------------------------------------------------------------------------- #
#                       Check if is MRSA                                    #

	nucl_db_dir = create_dir(raw_data, "nucl_db_dir")
	nucl_db = makeblastdb(makeblastdb_exe, nucl_db_dir, ffn, "nucl", "nucl_db")
	nucl_db_path = os.path.join(nucl_db_dir, nucl_db)
	prot_db_dir = create_dir(raw_data, "prot_db_dir")
	prot_db = makeblastdb(makeblastdb_exe, prot_db_dir, faa, "prot", "prot_db")
	prot_db_path = os.path.join(prot_db_dir, prot_db)

	print("-"*78)
	print("-"*78)
	print("-"*78)
	print

	orfx_nucl_hit = blast_parser(blastn_exe, nucl_db_path, orfx_base, "orfX")
	print("orfX Hit: ", orfx_nucl_hit)
	print
	mecA_hit = blast_parser(blastp_exe, prot_db_path, mecA_base, "mec")
	print("mec Hit: ", mecA_hit)
	print
	ccr_hit = blast_parser(blastp_exe, prot_db_path, ccr_base, "ccr")
	print("ccr Hit: ", ccr_hit)

	print
	print("-"*78)
	print

	core_elements = [orfx_nucl_hit, mecA_hit, ccr_hit]

#

	if all(core_elements):
		print("Buscando Contig ID de orfX, mec y ccr: ")
		print
		contig_id_orfx = get_contig(gff, orfx_nucl_hit)
		contig_id_mecA = get_contig(gff, mecA_hit)
		contig_id_ccr = get_contig(gff, ccr_hit)
		contig_ids = [contig_id_orfx, contig_id_mecA, contig_id_ccr]
		print("orfx in: ", contig_id_orfx)
		print("mec in: ", contig_id_mecA)
		print("ccr in: ", contig_id_ccr)
		print('\n'+"-"*78+'\n')
	else:
		print("Not MRSA")
		sys.exit()


	contigs_dict = fasta2dict(fna)
	nucl_dict = fasta2dict(ffn)
	actual_orfx = get_sequence(nucl_dict, orfx_nucl_hit)
	actual_mecA = get_sequence(nucl_dict, mecA_hit)
	actual_ccr = get_sequence(nucl_dict, ccr_hit)
	

	print actual_orfx

	if checkContig(contig_ids):

		# ------------------------------------------------------------------------- #
		#    Check if sccmec core components are in the same contig to continue     #
		# Template DNA contiene los genes orfx, mec y ccr (se puede usar contig_id_ orfx, ccr o mecA)

		template_dna = get_sequence(contigs_dict, contig_id_orfx)
		print("Largo TEMPLATE DNA: ", len(template_dna))

		# ------------------------------------------------------------------------- #
		# Invertir o Conservar sentido de la secuencia para que attL quede a la izq #

		template_dna = checkSense(actual_orfx, template_dna)


		texto_orfx, inicio_orfx, final_orfx = find_position(template_dna, actual_orfx, "orfX")

		print("ORFX INFO: ", texto_orfx, inicio_orfx, final_orfx)
		att_actual_orfx = actual_orfx[len(actual_orfx)-20:]
		print("ATT ACTUAL ORFX: ", att_actual_orfx)



		attL_database_dir = create_dir(raw_data, "attL_template_db")

		print("CURRENT DIR: ", os.getcwd())
		os.chdir(attL_database_dir)
		print("CURRENT DIR: ", os.getcwd())

		attL_file = "attL_{0}.fasta".format(nombre)
		print("ATTL FILE NAME: ", attL_file)
		with open(attL_file, "w") as f:
			f.write(">attL_"+nombre+"\n")
			f.write(str(att_actual_orfx+"\n"))


		attL_path = os.path.join(attL_database_dir, attL_file)
		database = makeblastdb(makeblastdb_exe, attL_database_dir, attL_path, "nucl", "attL_nucl_db")



		sccmec_right = template_dna[final_orfx:]
		with open("sccmec_RIGHT_"+nombre+".fasta", "w") as f:
			f.write(">sccmec_"+nombre+"\n")
			f.write(sccmec_right+"\n")
		print len(sccmec_right)

		size_window = len(att_actual_orfx)

		windows = ["".join(x) for x in window(sccmec_right, size_window)]

		core = "TATCATAA"


		# for i in windows:
		# 	if core in i:
		# 		result, err = blast(blastn_exe, database, i)
		# 		#print result
		# 		for hit in result.split('\n'):
		# 			if hit:
		# 				#print hit
		# 				args = [arg for arg in hit.split('\t')]
		# 				porcentaje = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100
		# 				if porcentaje >= 60.0:
		# 					print hit, porcentaje
		# 					print i 
		# 					m = re.search(i, sccmec_right)
		# 					n = re.search(reverse_complement(i), sccmec_right)
		# 					print("Same Sense: ", m) 
		# 					print("Anti Sense: ", n)
		# 					print("\n"+"-"*12+"\n")


		if re.findall(core, sccmec_right):
			found = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(core), sccmec_right)]
		print('\n')
		for f in found:
			print sccmec_right[f[0]-8:f[1]+4]
			print att_actual_orfx
			print('\n')

		# Usage
		print get_score('murat', 'bukar')

if __name__ == '__main__':
        main()