import sys
import os
import subprocess
import pandas as pd

def prokka_files_sccmec(location):
	""" ubicacion de archivos prokka a utilizar """
	for file in os.listdir(location):
		if file.endswith(".ffn"):
			ffn = os.path.join(location, file)
		if file.endswith(".gff"):
			gff = os.path.join(location, file)
		if file.endswith(".fna"):
			fna = os.path.join(location, file)
		if file.endswith(".faa"):
			faa = os.path.join(location, file)
		if file.endswith(".gbk"):
			gbk= os.path.join(location, file)

	return ffn, gff, fna, faa, gbk


# ------------------------------------------------------------------------- #

def execute_prokka_sccmec(prokka_exe, output_prokka, contigs):
	""" RUN PROKKA parametros especiales para s aureus """
	#print output_prokka
	kingdom = "Bacteria"
	genus = "Staphylococcus"
	locustag = output_prokka
	cmd_prokka = prokka_exe+" --kingdom " + kingdom + \
					" --outdir " + output_prokka + \
					" --quiet "+ \
					" --genus "+ genus + \
					" --locustag "+ locustag + \
					" --centre 10 --compliant" + \
					" --prefix "+ output_prokka + \
					" --force " + \
					' '+ contigs
	os.system(cmd_prokka)


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

def sBLAST(blast_exe, database, query):
	""" Ejecuta BLAST usando database+string """
	outfmt = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"

	cmd = [blast_exe, "-outfmt", outfmt, "-db", database]

	proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
								stdout=subprocess.PIPE,
								stderr=subprocess.STDOUT)
	
	results, err = proc.communicate(query)

	return results, err

# ------------------------------------------------------------------------- #

def annotation_blast_parser(out):
	""" Get Best Hit If Any """
	salida = ""
	if out:
		row = [s.split('\t') for s in out.split('\n')]
		lines = [x for x in row if x != [""]]
		for line in lines:
			args = [arg for arg in line]

			porcentaje_query = float("{0:.2f}".format(((float(args[5])-float(args[4]))+1)/float(args[1])*100))
			porcentaje_subject = float("{0:.2f}".format(((float(args[7])-float(args[6]))+1)/float(args[3])*100))

			#print args, porcentaje_query, porcentaje_subject

			if porcentaje_query >= 70.0 and porcentaje_subject >= 40.0:

				args.append(str(porcentaje_query))
				args.append(str(porcentaje_subject))
				salida += "\t".join([x for x in args])+'\n'
	else:
		return None

	if salida:
		lines = [s.split('\t') for s in salida.split('\n')]
		list2 = [x for x in lines if x != [""]]

		list2.sort(key=lambda x: float(x[10]), reverse=True)
		
		#print("Best Hit: ", list2[0])
		if float(list2[0][13]) < 60.0:
			saureus_id = list2[0][2]+"-Truncated"
			return saureus_id
		else:

			saureus_id = list2[0][2]
			return saureus_id 

	else:

		return None

# ------------------------------------------------------------------------- #

def locationGenbank(genbank):
    db_dict_pos = {}
    with open(genbank) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            if " gene " in lines[i]:
                #print len(lines[i])
                #print lines[i].strip()
                posicion = lines[i].strip()
                posicion = posicion[16:]
                
            if "/locus_tag" in lines[i]:
                #print lines[i].strip()
                locus = lines[i].strip().split('"')[1]
                db_dict_pos[posicion] = locus


    ordered_keys = []
    for key in db_dict_pos.keys():
        ordered_keys.append((db_dict_pos[key], key))

    ordered_keys = sorted(ordered_keys)
    new_dict = {}
    for i in ordered_keys:
        locus = i[0]
        if "complement" in i[1]:
            #print i[1].replace(")", "").replace("complement(", "(-)\t").replace("..", "\t")
            loc = i[1].replace(")", "").replace("complement(", "(-)\t").replace("..", "\t")
        else:
            #print "(+)\t"+i[1].replace("..", "\t")
            loc = "(+)\t"+i[1].replace("..", "\t")
        
        sense, start, end = loc.split()
        largo_nucl = int(end)-int(start)
        
        #print locus, sense, start, end, largo_nucl
        
        info = sense+'\t'+start+'\t'+end+'\t'+str(largo_nucl)
        #print info
        new_dict[locus] = info

    return new_dict

# ------------------------------------------------------------------------- #

ccrA1B1 = ["ccrA1", "ccrB1"]
ccrA1B1_ = ["ccrA1", "ccrB1-Truncated"]
ccrA1B1__ = ["ccrA1-Truncated", "ccrB1"]

ccrA2B2 = ["ccrA2", "ccrB2"]
ccrA2B2_ = ["ccrA2", "ccrB2-Truncated"]
ccrA2B2__ = ["ccrA2-Truncated", "ccrB2"]

ccrA3B3 = ["ccrA3", "ccrB3"]
ccrA3B3_ = ["ccrA3", "ccrB3-Truncated"]
ccrA3B3__ = ["ccrA3-Truncated", "ccrB3"]

ccrA4B4 = ["ccrA4", "ccrB4"]
ccrA4B4_ = ["ccrA4", "ccrB4-Truncated"]
ccrA4B4__ = ["ccrA4-Truncated", "ccrB4"]

ccrC = ["ccrC"]

ccrA1B3 = ["ccrA1", "ccrB3"]

def ccrAllotype(lst):
	if all(gene in lst for gene in ccrA1B1):
		return "1"
	if all(gene in lst for gene in ccrA1B1_):
		return "1"	
	if all(gene in lst for gene in ccrA1B1__):
		return "1"

	if all(gene in lst for gene in ccrA2B2):
		return "2"
	if all(gene in lst for gene in ccrA2B2_):
		return "2"
	if all(gene in lst for gene in ccrA2B2__):
		return "2"

	if all(gene in lst for gene in ccrA3B3):
		return "3"
	if all(gene in lst for gene in ccrA3B3_):
		return "3"
	if all(gene in lst for gene in ccrA3B3__):
		return "3"

	if all(gene in lst for gene in ccrA4B4):
		return "4"
	if all(gene in lst for gene in ccrA4B4_):
		return "4"
	if all(gene in lst for gene in ccrA4B4__):
		return "4"

	if all(gene in lst for gene in ccrC):
		return "5"	

	if all(gene in lst for gene in ccrA1B3):
		return "8"

	else:
		return "N"


def mecClass2(lst):
	for i in range(len(lst)):
		if lst[i] == "mecA":

			if lst[i+1] == "mecR1":
				if lst[i+2] == "MecI":
					return "A"
			if lst[i+1] == "mecR1":
				if lst[i+2] == "MecI-Truncated":
					return "A"
			if lst[i+1] == "mecR1-Truncated":
				if lst[i+2] == "MecI":
					return "A"
			if lst[i+1] == "mecR1-Truncated":
				if lst[i+2] == "MecI-Truncated":
					return "A"

			if lst[i+1] == "mecR1-Truncated":
				if lst[i+2] == "IS1272":
					if lst[i-1] == "IS431":
						return "B"


			if lst[i-1] == "IS431":
				if lst[i+1] == "mecR1-Truncated":
					if lst[i+2] == "IS431":
						return "C1"

			if lst[i-1] == "IS431" or lst[i-1] == "IS431-Truncated":
				if lst[i+1] == "IS431" or lst[i+1] == "IS431-Truncated":

					try:
						if "ccrA" in lst[i-2] or "ccrB" in lst[i-2]:
							return "CL"
						if "ccrA" in lst[i+2] or "ccrB" in lst[i+2]:
							return "CR"
					except IndexError:
						pass

				return "C"


		if lst[i] == "mecC":
			if lst[i+1] == "mecR1":
				if lst[i+2] == "MecI":
					return "E"

	else:
		return "N"


sccmec_dict = {
	"1B": "I",
	"2A": "II",
	"3A": "III",
	"2B": "IV",
	"5C": "V",
	"4B": "VI",
	"5C1": "VII",
	"4A": "VIII",
	"1CL": "IX",
	"1CR": "X",
	"8E": "XI"
}
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #


def current_annotation(sccmec_file, blast_exe, prokka_exe, core_database):

	sccmec_id = sccmec_file.split(".")[0].split('/')[-1]

	execute_prokka_sccmec(prokka_exe, sccmec_id, sccmec_file)
	ffn, gff, fna, faa, gbk = prokka_files_sccmec(sccmec_id)
	faa_dict = fasta2dict(faa)

	position_dict = locationGenbank(gbk)

	print("\n"+"-"*42)
	print sccmec_id

	ordered_keys = []
	for key in faa_dict.keys():
		ordered_keys.append(key)


	ordered_keys = sorted(ordered_keys)

	locus, sentido, inicio, fin, nucl, prot, gene = ([] for i in range(7))
	for k in ordered_keys:
		sequence = get_sequence(faa_dict, k)
		largo = str(len(sequence))
		results, err = sBLAST(blast_exe, core_database, sequence)
		match = annotation_blast_parser(results)

		#print position_dict[k]
		sense, start, end, size = position_dict[k].split()

		if match:
			match = match.split("_")[-1]
			#print k, sense, start, end, size, largo, match
			locus.append(k)
			sentido.append(sense)
			inicio.append(start)
			fin.append(end)
			nucl.append(size)
			prot.append(largo)
			gene.append(match)

		else:
			match = "NN"
			#print k, sense, start, end, size, largo, match
			locus.append(k)
			sentido.append(sense)
			inicio.append(start)
			fin.append(end)
			nucl.append(size)
			prot.append(largo)
			gene.append(match)
			continue

	columns = ["Locus", "Sense", "Start", "End", "Size", "Length", "Gene"]
	data = pd.DataFrame({
		"Locus": locus,
		"Sense": sentido,
		"Start": inicio,
		"End": fin,
		"Size": nucl,
		"Length": prot,
		"Gene": gene
		}, columns=columns)

	os.chdir(sccmec_id)

	file_name = "annotation_"+sccmec_id+".txt"
	data.to_csv(file_name, sep='\t', encoding='utf-8', index=False)

	gene = [x for x in gene if x != "NN"]

	print gene
	core_elements = pd.DataFrame({"Core Elements": gene})
	file_name = "core_elements_"+sccmec_id+".txt"
	core_elements.to_csv(file_name, sep='\t', encoding='utf-8', index=False)

	clase = mecClass2(gene)
	alotipo = ccrAllotype(gene)
	code = (alotipo+clase)

	try: 
		print sccmec_dict[code]
		file_name = "type_"+sccmec_id+".txt"
		with open(file_name, "w") as f:
			f.write(sccmec_dict[code]+'\n')
	except KeyError:
		print code
		file_name = "type_"+sccmec_id+".txt"
		with open(file_name, "w") as f:
			f.write(code+'\n')

