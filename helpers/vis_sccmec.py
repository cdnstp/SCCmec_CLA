import sys
import subprocess
from dna_features_viewer import (GraphicFeature, GraphicRecord)

colors = {
	'orfX': 'green', 'IS431': 'blue', 'IS1272': 'blue',
	'mecA': 'red', 'mecR1': 'red', 'MecI': 'red',
	'mecB': 'red', 'mecC': 'red', 'ccrC': 'yellow', 
	'ccrA1': 'yellow', 'ccrA2': 'yellow', 'ccrA2': 'yellow', 
	'ccrA3': 'yellow', 'ccrA4': 'yellow', 'ccrA7': 'yellow', 
	'ccrB1': 'yellow', 'ccrB1': 'yellow', 'ccrB2': 'yellow', 
	'ccrB2': 'yellow', 'ccrB3': 'yellow', 'ccrB4': 'yellow', 
	'NN': '#0d0d0d', 'mecR1-Truncated': 'red', 'MecI-Truncated': 'red',
}

labels = {
	'NN': None,
	'mecR1-Truncated': 'mecR1-Truncated',
	'MecI-Truncated': 'mecI-Truncated',
	'MecI': 'mecI',
	'mecR1': 'mecR1',
	'ccrB3': 'ccrB3',
	'ccrB2': 'ccrB2',
	'ccrB1': 'ccrB1',
	'mecC': 'mecC',
	'mecB': 'mecB',
	'mecA': 'mecA',
	'orfX': 'orfX',
	'ccrB4': 'ccrB4',
	'IS1272': 'IS1272',
	'IS431': 'IS431',
	'ccrC': 'ccrC',
	'ccrA2': 'ccrA2',
	'ccrA3': 'ccrA3',
	'ccrA1': 'ccrA1',
	'ccrA7': 'ccrA7',
	'ccrA4': 'ccrA4',
}

# ------------------------------------------------------------------------- #
def get_sequence(sequences_dict, sequence_id):
	""" Obtiene secuencia desde dict-seq usando la key """
	sequence = ""
	for i in sequences_dict.get(sequence_id)[:]:
		sequence += ''.join(i)
	return sequence

def fasta2dict(file):
	""" Transforma archivos multi-fasta en un dict de python """
	fastadict = {}
	with open(file) as file_one:
		for line in file_one:
			line = line.strip()
			if not line:
				continue
			if line.startswith('>'):
				active_sequence_name = str(line.split()[0][1:])
				if active_sequence_name not in fastadict:
					fastadict[active_sequence_name] = []
				continue
			sequence = line
			fastadict[active_sequence_name].append(sequence)
	return fastadict

def run_blastp(blastp, subject, query):
	formato = '6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue'
	process = subprocess.Popen([blastp, '-subject', subject, '-outfmt', formato], 
				stdin=subprocess.PIPE,
				stdout=subprocess.PIPE,
				stderr=subprocess.STDOUT)
	out, err = process.communicate(query)

	return out, err

# ------------------------------------------------------------------------- #
# guardar annotation file en una lista para evaluar si genes anotados como
# NN estan presentes en core-proteins de los cassettes de este grupo
# y asi cambiar la visualizacion
def annotation_data(file):
	datafile = []
	with open(file) as f:
		lines = f.readlines()
		for line in lines:
			line = line.strip()
			if not line: continue
			if 'Locus' in line: continue
			datafile.append(line)
	return datafile

# ------------------------------------------------------------------------- #
def output_blastp(out):
	besthit = []
	hits = [hit.split('\t') for hit in out.split('\n')]
	for hit in hits:
		if [''] == hit: continue
		coverage = (((float(hit[5])-float(hit[4]))+1)/float(hit[1]))*100
		if coverage >= 70.0 and float(hit[10]) >= 50.0:
			besthit.append(hit)

	if besthit:
		return besthit[0]
	else:
		return None


def update_annotation(lista, blastp, faa_dict_sccmec, core_proteins):
	update_datafile = []
	for line in lista:
		print(line)
		id_, sense, start, end, size, length, gene = line.split()

		if 'NN' == gene:
			seq = get_sequence(faa_dict_sccmec, id_)
			fastaformat = '>{0}\n{1}\n'.format(id_, seq)
			out, err = run_blastp(blastp, core_proteins, fastaformat)
			if out:
				besthit = output_blastp(out)
				if besthit:
					gene = 'core-proteins'
					update_datafile.append([id_, sense, start, end, size, length, gene])
				else:
					update_datafile.append([id_, sense, start, end, size, length, gene])
			else:
				update_datafile.append([id_, sense, start, end, size, length, gene])
		else:
			update_datafile.append([id_, sense, start, end, size, length, gene])

	return update_datafile


# ------------------------------------------------------------------------- #

def vis_sccmec(faa_file_sccmec, annotation_file, length_sccmec, core_proteins, blastp):

	# use faa file from prokka annotation on sccmec
	faa_dict_sccmec = fasta2dict(faa_file_sccmec)

	# update annotation based on core proteins in cluster
	datafile = annotation_data(annotation_file)
	update_datafile = update_annotation(datafile, blastp, faa_dict_sccmec, core_proteins)

	# create features object to visualisation using dna_features_viewer
	features = []
	for line in update_datafile:
		id_, sense, start, end, size, length, gene = line
		if gene == 'core-proteins':
			color = '#ff8848'
			label = None
		else:
			try:
				color = colors[gene]
			except KeyError:
				color = 'grey'
			try:
				label = labels[gene]
			except KeyError:
				label = None
		if '-' in sense:
			features.append(GraphicFeature(start=int(start), end=int(end), strand=-1,
											color=color, label=label))

		if '+' in sense:
			features.append(GraphicFeature(start=int(start), end=int(end), strand=+1,
											color=color, label=label))

	record = GraphicRecord(sequence_length=length_sccmec, features=features)
	ax, _ = record.plot(figure_width=20)

	id_ = annotation_file.split('_')[-1].split('.')[0]
	filename = 'SCCmec_{}.png'.format(id_)

	ax.figure.savefig(filename, dpi=300)








