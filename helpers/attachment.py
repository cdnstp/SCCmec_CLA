from blast_helpers import *
import re

# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #

def get_attchment(nombre, output, template_dna, attr_database_path, blastn_exe):
	os.chdir(output)
	core = "TATCATAA"
	print("BUSCANDO ATTR")
	matches = re.finditer("{}".format(core), template_dna)
	salida = ""
	atts_location = {}
	for match in matches:
		span = match.span()
		desde = span[0]-28
		hasta = span[1]+22

		# Para evitar considerar attL2, posicion de attr encontrado > 1000pb
		if hasta >= 1000 and hasta <= 100000:
			query = template_dna[desde:hasta]
			atts_location[query] = [desde, hasta]
			out, err = attR_BLAST(blastn_exe, attr_database_path, query)

			row = [s.split('\t') for s in out.split('\n')]
			lines = [x for x in row if x != [""]]

			for line in lines:

				args = [arg for arg in line]
				porcentaje = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100
				#print line, porcentaje
				if porcentaje >= 30.0:
					args.append(str(format(porcentaje, '.2f')))
					args.append(str(desde))
					args.append(str(hasta))
					args[0] = query
					salida += "\t".join([x for x in args])+'\n'

	#print salida
	#print atts_location
	#print type(atts_location)
	#print len(atts_location)

	if len(atts_location) == 1:
		print("1 attR encontrado")
		attR_sequence = atts_location.keys()[0]
		attR_inicio = atts_location.values()[0][0]
		attR_final = atts_location.values()[0][1]
		hit = "attR_"+nombre+"_"+str(attR_inicio)+"-"+str(attR_final)
		with open("attR_"+nombre+".fasta", "w") as f:
			f.write(">"+hit+"\n")
			f.write(attR_sequence+"\n")

		return int(attR_inicio), int(attR_final), hit

	if salida:
		lines = [s.split('\t') for s in salida.split('\n')]
		list2 = [x for x in lines if x != [""]]
		list2.sort(key=lambda x: float(x[10]), reverse=True)
		print("Best Hit attR: ", list2[0])
		attR_sequence = list2[0][0]

		print attR_sequence
		attR_inicio = list2[0][13]
		attR_final = list2[0][14]
		hit = "attR_"+nombre+"_"+str(attR_inicio)+"-"+str(attR_final)
		with open("attR_"+nombre+".fasta", "w") as f:
			f.write(">"+hit+"\n")
			f.write(attR_sequence+"\n")

		return int(attR_inicio), int(attR_final), hit

	else:
		print("Not even TATCATAA was found")
		return None
