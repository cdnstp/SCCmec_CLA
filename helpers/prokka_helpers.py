import os

# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #

def prokka_files(location):
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

	return ffn, gff, fna, faa


# ------------------------------------------------------------------------- #

def execute_prokka(prokka_exe, output_prokka, contigs, name):
	""" RUN PROKKA parametros especiales para s aureus """
	#					" --force " + \
	# falta agregar --force
	kingdom = "Bacteria"
	genus = "Staphylococcus"

	prefix_locus = name
	cmd_prokka = prokka_exe+" --kingdom " + kingdom + \
					" --outdir " + output_prokka + \
					" --quiet "+ \
					" --genus "+ genus + \
					" --locustag "+ prefix_locus + \
					" --centre 10 --compliant" + \
					" --prefix "+ prefix_locus + \
					' '+ contigs
	os.system(cmd_prokka)


# ------------------------------------------------------------------------- #

def get_region(gff, query):
	""" Entrega el ID del contig en el que esta el query """
	with open(gff) as f:
		header = []
		lines = f.readlines()
		for line in lines:
			if query in line:
				header.append(line.split()[0])

	return header[0]


# ------------------------------------------------------------------------- #

