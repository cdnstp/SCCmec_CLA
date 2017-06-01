import os
import subprocess

# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #

def makeblastdb(makeblastdb_exe, working_dir, input_seq, nucl, db_name):
	os.chdir(working_dir)
	proc = makeblastdb_exe+' -in '+input_seq+' -dbtype '+nucl+' -out '+db_name
	os.system(proc)

	return db_name


# ------------------------------------------------------------------------- #

def simpleBlast(blast_exe, database, query, qname):
	""" Ejecuta BLAST usando database + file """
	formato = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"

	process = subprocess.Popen([blast_exe, "-query", query, 
								"-db", database, 
								"-outfmt", formato], stdin=subprocess.PIPE,
													stdout=subprocess.PIPE,
													stderr=subprocess.STDOUT)

	out, err = process.communicate()

	return out, err


# ------------------------------------------------------------------------- #

def simpleBlastParser(out):
	""" Get Best Hit If Any """
	salida = ""
	if out:
		row = [s.split('\t') for s in out.split('\n')]
		lines = [x for x in row if x != [""]]
		for line in lines:
			args = [arg for arg in line]
			porcentaje = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100
			if porcentaje >= 70.0:
				args.append(str(porcentaje))
				salida += "\t".join([x for x in args])+'\n'
	else:
		return None

	if salida:
		lines = [s.split('\t') for s in salida.split('\n')]
		list2 = [x for x in lines if x != [""]]

		list2.sort(key=lambda x: float(x[10]), reverse=True)
		print("Best Hit: ", list2[0])

		saureus_id = list2[0][2]

		if float(list2[0][10]) >= 50.0:
			return saureus_id 

		else:
			return None

	else:
		return None


# ------------------------------------------------------------------------- #

def attR_BLAST(blast_exe, database, query):
	""" Ejecuta BLAST usando database+string """
	outfmt = "6 qseqid qlen sseqid slen qstart qend sstart send length nident pident evalue"

	cmd = [blast_exe, "-outfmt", outfmt, "-word_size", "6", "-db", database]

	proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
								stdout=subprocess.PIPE,
								stderr=subprocess.STDOUT)
	
	results, err = proc.communicate(query)

	return results, err


# ------------------------------------------------------------------------- #