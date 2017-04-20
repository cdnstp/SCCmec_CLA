def findAtts(nombre, output, sequence, att):
	os.chdir(output)
	"""?att patron basico corresponde a los nucl finales del gen orfX"""
	atts_location = {}
	for i in re.findall("({s}){{s<=5}}".format(s=att), sequence):
		atts_location[seq[m.start(0):(m.end(0))]] = [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), seq)]

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

	for s in sites:
		fr = s[1][0][0]
		to = s[1][0][1]
		if fr+to > (val_attl+300):
			hipotetical_attR[s[0]] = s[1][0]
		else:
			hipotetical_attL[s[0]] = s[1][0]

		attL_list = sorted(hipotetical_attL.items(), key=lambda e: e[1][0])
		attR_list = sorted(hipotetical_attR.items(), key=lambda e: e[1][0])

		attL_sequence = seq[attL_list[0][1][0]-20:attL_list[0][1][1]+20]
		with open("attL_"+nombre+".fasta", "w") as f:
			f.write(">attL_"+nombre+"\n")
			f.write(attL_sequence+'\n')

		with open("attR_"+nombre+".fasta", "w") as f:
			for attr in attR_list:
				#print attr
				f.write(">attR_"+nombre+"_"+str(attr[1][0]-20)+"_"+str(attr[1][1]+20)+"\n")
				f.write(str(seq[attr[1][0]-20:attr[1][1]+20]))
				#print seq[attr[1][0]-20:attr[1][1]+20]

	coordinates = [attL_list[0][1][0], attr[-1][1][1]]

	return coordinates

import subprocess
import os
from subprocess import Popen, PIPE

blast_exe = "/usr/bin/blastn"
attr_db = "/home/fsepulveda/Desktop/test_script/input/attR/attR_database.txt"
attr = "/home/fsepulveda/Desktop/test_script/output_CP007676/attR_CP007676.fasta"

def blastAlign(blast_exe, query, subject):
	formato = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"
	process = subprocess.Popen([blast_exe, "-query", query, "-subject", subject, "-outfmt", formato], stdin=subprocess.PIPE,
		stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT)
	out, err = process.communicate()

	return out, err


result, err = blastAlign(blast_exe, attr, attr_db)

hit = result.split()[0]
print hit
from_loc = hit.split("_")[-1].split("-")[0]
to_loc = hit.split("_")[-1].split("-")[1]

print "Desde posicion: ", from_loc, "hasta: ", to_loc
