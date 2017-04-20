import os
import sys
import regex as re
import subprocess




def config():
	import ConfigParser
	config = ConfigParser.ConfigParser()
	config.readfp(open(r'config.txt'))
	prokka_exe = config.get('configuration', 'prokka')
	blastn_exe = config.get('configuration', 'blastn')
	blastp_exe = config.get('configuration', 'blastp')
	makeblastdb_exe = config.get('configuration', 'makeblastdb')
	inputFiles = config.get('configuration', 'input')
	contigs = config.get('configuration', 'contig')
	ccr = config.get('configuration', 'ccr')
	orfX = config.get('configuration', 'orfX')
	mecA = config.get('configuration', 'mecA')
	attr_db = config.get('configuration', 'attr_db')

	return prokka_exe, blastn_exe, blastp_exe, makeblastdb_exe, inputFiles, contigs, ccr, orfX, mecA, attr_db


def simple_sequence(file):
    with open(file) as f:
        lines = f.read()
        sequences = [line for line in lines.split('\n') if not ">" in line]
    sequence = ''.join(map(str,sequences))

    return sequence

def blastAlign(blast_exe, query, subject):
	formato = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"
	process = subprocess.Popen([blast_exe, "-query", query, "-subject", subject, "-outfmt", formato], stdin=subprocess.PIPE,
		stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT)
	out, err = process.communicate()
	
	return out, err

def blast(blast_exe, database, query):
	outfmt = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"
	cmd = [blast_exe, "-outfmt", outfmt, "-db", database]
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


def blast_parser(blast_exe, db, query):
	blastfile, err = blast(blast_exe, db, query)
	best_hit = blastfile.split()[2]

	return best_hit


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
	if re.findall("{pattern}".format(pattern=sequence), contig):
		return contig, "+"
	contig = reverse_complement(contig)
	if re.findall("{pattern}".format(pattern=sequence), contig):
		return contig, "+"


def create_dir(base, dir_name):
	new_dir = os.path.join(base, dir_name)
	if not os.path.isdir(new_dir):
		os.makedirs(new_dir)
	return new_dir


def find_position(sequence, gene, name):
	if re.findall("{pattern}".format(pattern=gene), sequence):
		for i in re.findall("{pattern}".format(pattern=gene), sequence):
			print "+ %s located at: " % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]
	else:
		for i in re.findall("{pattern}".format(pattern=reverse_complement(gene)), sequence):
			print "- %s located at: " % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]


def findAtts(nombre, output, sequence, att, attr_database_path, blastn_exe):
	os.chdir(output)
	"""?att patron basico corresponde a los nucl finales del gen orfX"""
	atts_location = {}
	for i in re.findall("({s}){{s<=5}}".format(s=att), sequence):
		atts_location[sequence[m.start(0):(m.end(0))]] = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]

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
	attL_pos = int([x+y for (x,y) in sites[0][1]][0])

	hipotetical_attL = {}
	hipotetical_attR = {}

	val_attl = int([x+y for (x,y) in sites[0][1]][0])

	for s in sites:
		fr = s[1][0][0]
		to = s[1][0][1]
		if fr+to > (val_attl+300):
			hipotetical_attR[s[0]] = s[1][0]
		else:
			hipotetical_attL[s[0]] = s[1][0]

		attL_list = sorted(hipotetical_attL.items(), key=lambda e: e[1][0])
		attR_list = sorted(hipotetical_attR.items(), key=lambda e: e[1][0])

		attL_sequence = sequence[attL_list[0][1][0]-20:attL_list[0][1][1]+20]
		with open("attL_"+nombre+".fasta", "w") as f:
			f.write(">attL_"+nombre+"\n")
			f.write(attL_sequence+'\n')
		output_attR = "attR_%s.fasta" % nombre
		with open(output_attR, "w") as f:
			for attr in attR_list:
				#print attr
				f.write(">attR_"+nombre+"_"+str(attr[1][0]-20)+"-"+str(attr[1][1]+20)+"\n")
				f.write(str(sequence[attr[1][0]-20:attr[1][1]+20])+'\n')
				#print seq[attr[1][0]-20:attr[1][1]+20]

	"""ADDING CHECK UP FOR ATTR """
	print os.getcwd()
	attR_list_path = os.path.join(output, output_attR)
	print attR_list_path
	result, err = blastAlign(blastn_exe, attR_list_path, attr_database_path)
	hit = result.split()[0]
	if hit:
		attR_start = hit.split("_")[-1].split("-")[0]
		attR_end = hit.split("_")[-1].split("-")[1]
		#print "attL starting at: ", attL_list[0][1][0]
		#print "attR ends at: ", attR_list[-1][1][1]
		coordinates = [attL_list[0][1][0], int(attR_end)]
		print coordinates

		return coordinates, hit


# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #


def main():
	prokka_exe, blastn_exe, blastp_exe, makeblastdb_exe, inputFiles, contigs, ccr, orfX, mecA, attr_db = config()

	orfx_base = simple_sequence(os.path.join(inputFiles, orfX))
	mecA_base = simple_sequence(os.path.join(inputFiles, mecA))
	ccr_base = simple_sequence(os.path.join(inputFiles, ccr))
	attr_database_path = os.path.join(inputFiles, attr_db)

	contigs = sys.argv[1]
	nombre = contigs.split("/")[-2].split("_")[1]
	
	#print nombre
	contigs = os.path.join(inputFiles, contigs)
	working_dir = os.getcwd()
	#print working_dir
	#output = "output_" + id_generator()
	output = "output_"+nombre
	raw_data = os.path.join(working_dir, output)
	output_prokka = os.path.join(raw_data, "output_prokka")

	# execute prokka 
	execute_prokka(prokka_exe, output_prokka, contigs)
	ffn, gff, fna, faa = prokka_files(output_prokka)

	# create nucleotide db
	nucl_db_dir = create_dir(raw_data, "nucl_db_dir")
	nucl_db = makeblastdb(makeblastdb_exe, nucl_db_dir, ffn, "nucl", "nucl_db")
	nucl_db_path = os.path.join(nucl_db_dir, nucl_db)
	# run blastn and parse output
	orfx_nucl_hit = blast_parser(blastn_exe, nucl_db_path, orfx_base)


	prot_db_dir = create_dir(raw_data, "prot_db_dir")
	prot_db = makeblastdb(makeblastdb_exe, prot_db_dir, faa, "prot", "prot_db")
	prot_db_path = os.path.join(prot_db_dir, prot_db)
	# run blastp for mecA and ccrs
	mecA_hit = blast_parser(blastp_exe, prot_db_path, mecA_base)
	ccr_hit = blast_parser(blastp_exe, prot_db_path, ccr_base)

	# find contig which contains orfx, mec and ccr
	contig_id_orfx = get_contig(gff, orfx_nucl_hit)
	contig_id_mecA = get_contig(gff, mecA_hit)
	contig_id_ccr = get_contig(gff, ccr_hit)
	contig_ids = [contig_id_orfx, contig_id_mecA, contig_id_ccr]


	print
	print("-"*78)
	print

	print contig_ids
	""" Check if sccmec core components are in the same contig to continue """
	if checkContig(contig_ids):

		contigs_dict = fasta2dict(fna)
		seq = get_sequence(contigs_dict, contig_id_orfx)
		nucl_dict = fasta2dict(ffn)
		actual_orfx = get_sequence(nucl_dict, orfx_nucl_hit)
		actual_mecA = get_sequence(nucl_dict, mecA_hit)
		actual_ccr = get_sequence(nucl_dict, ccr_hit)


		""" Check if orfX has the same sense as the contig or use the reverse complementary sequence
			of the given contig
		"""
		seq, orfx_sense = checkSense(actual_orfx, seq)

		""" Extract 19 nucleotides located at the 3'- end of orfX gene corresponding to attL """
		att_actual_orfx = actual_orfx[len(actual_orfx)-20:-2]
		print att_actual_orfx

		""" Search for attachment site sequences """
		""" Filter att sequences according to literature """

		coordinates, hit = findAtts(nombre, raw_data, seq, att_actual_orfx, attr_database_path, blastn_exe)
		sccmec = seq[coordinates[0]:coordinates[-1]]

		print "Contig Length: ", len(seq)
		print "SCCmec Length: ", len(sccmec)
		print "attR match: ", hit


		""" EDITAR """
		""" ADD POSITION CHECK """
		#find_position(seq, actual_mecA, "mecA")
		#find_position(seq, actual_ccr, "ccr")

		with open("sccmec_"+nombre+".fasta", "w") as f:
			f.write(">sccmec_"+nombre+"_"+str(len(sccmec))+"\n")
			for i in range(0, len(sccmec), 60):
				f.write(sccmec[i:i+60]+'\n')

	else:
		print "Error: core components are not in the same contig"
		sys.exit()


if __name__ == '__main__':
        main()