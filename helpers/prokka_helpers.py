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

def execute_prokka(prokka_exe, output_prokka, contigs):
	""" RUN PROKKA parametros especiales para s aureus """
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


# ------------------------------------------------------------------------- #

def get_contig(gff, query):
	""" Entrega el ID del contig en el que esta el query """
	with open(gff) as f:
		header = []
		lines = f.readlines()
		for line in lines:
			if query in line:
				header.append(line.split()[0])

	return header[0]


# ------------------------------------------------------------------------- #

