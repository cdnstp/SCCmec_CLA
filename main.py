from helpers import *


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


# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #


def main():

# ------------------------------------------------------------------------- #
#                      CONFIGURATION, OUTPUT SET UP                         #

	prokka_exe, blastn_exe, blastp_exe, makeblastdb_exe, inputFiles, seqFile, ccr, orfX, mecA, attr_db = config()

	orfx_base = os.path.join(inputFiles, orfX)
	mecA_base = os.path.join(inputFiles, mecA)
	ccr_base = os.path.join(inputFiles, ccr)

	attr_database_path = os.path.join(inputFiles, attr_db)

	nombre = seqFile.split(".")[0]
	print nombre
	contigs = os.path.join(inputFiles, seqFile)
	working_dir = os.getcwd()
	output = "output_"+nombre
	raw_data = os.path.join(working_dir, output)
	output_prokka = os.path.join(raw_data, "output_prokka")

	execute_prokka(prokka_exe, output_prokka, contigs)
	ffn, gff, fna, faa = prokka_files(output_prokka)

	contigs_dict = fasta2dict(fna)
	nucl_dict = fasta2dict(ffn)

# ------------------------------------------------------------------------- #
#                       Check if is MRSA                                    #

# ------------------------------------------------------------------------- #
#                       Database Setup                                      # 

	nucl_db_dir = create_dir(raw_data, "nucl_db_dir")
	nucl_db = makeblastdb(makeblastdb_exe, nucl_db_dir, ffn, "nucl", "nucl_db")
	nucl_db_path = os.path.join(nucl_db_dir, nucl_db)

	prot_db_dir = create_dir(raw_data, "prot_db_dir")
	prot_db = makeblastdb(makeblastdb_exe, prot_db_dir, faa, "prot", "prot_db")
	prot_db_path = os.path.join(prot_db_dir, prot_db)

	print("-"*78+"\n"+"-"*78+"\n"+"-"*78+"\n")

# ------------------------------------------------------------------------- #
#                 BLAST orfX + PBP2a + ccr                                  #

	print("\t"+"    RUNNING BLAST    "+"\n")

	print orfx_base
	orfx_nucl_hit, err = simpleBlast(blastn_exe, nucl_db_path, orfx_base, "orfX")
	orfx_nucl_hit = simpleBlastParser(orfx_nucl_hit)
	print("orfX Hit: ", orfx_nucl_hit)
	print("\n"+"-"*78+"\n")

	print mecA_base
	mecA_hit, err = simpleBlast(blastp_exe, prot_db_path, mecA_base, "PBP2a")
	mecA_hit = simpleBlastParser(mecA_hit)
	print("PBP2a Hit: ", mecA_hit)
	print("\n"+"-"*78+"\n")

	print ccr_base
	ccr_hit, err = simpleBlast(blastp_exe, prot_db_path, ccr_base, "ccr")
	ccr_hit = simpleBlastParser(ccr_hit)
	print("ccr Hit: ", ccr_hit)
	print("\n"+"-"*78+"\n")

	core_elements = [orfx_nucl_hit, mecA_hit, ccr_hit]


	if all(core_elements):

		print("CORE ELEMENTS ARE PRESENT "+"\n")


		contig_id_orfx = get_contig(gff, orfx_nucl_hit)
		contig_id_mecA = get_contig(gff, mecA_hit)
		contig_id_ccr = get_contig(gff, ccr_hit)

		contig_ids = [contig_id_orfx, contig_id_mecA, contig_id_ccr]

		print("Contig orfX: ", contig_id_orfx)
		print("Contig PBP2a: ", contig_id_mecA)
		print("Contig ccr: ", contig_id_ccr)
		print("\n"+"-"*78+"\n")


	else:
		print("Not MRSA")
		os.chdir(raw_data)
		with open("info_"+nombre+".txt", "w") as f:
			f.write("Not MRSA"+'\n')
			f.write("OrfX: "+str(core_elements[0])+'\n')
			f.write("PBP2a: "+str(core_elements[1])+'\n')
			f.write("CCR: "+str(core_elements[2])+'\n')
		if core_elements[0] is not None:
			actual_orfx = get_sequence(nucl_dict, orfx_nucl_hit)
			att_actual_orfx = actual_orfx[len(actual_orfx)-21:]
			with open("attB_"+nombre+".fasta", "w") as f:
				f.write(">attB_"+nombre+"\n")
				f.write(att_actual_orfx+"\n")	
		sys.exit()

# ------------------------------------------------------------------------- #
#            Check if CORE ELEMENTS are in the same contig                  #

	actual_orfx = get_sequence(nucl_dict, orfx_nucl_hit)
	actual_mecA = get_sequence(nucl_dict, mecA_hit)
	actual_ccr = get_sequence(nucl_dict, ccr_hit)
	
	if checkContig(contig_ids):
		os.chdir(raw_data)
		template_dna = get_sequence(contigs_dict, contig_id_orfx)

		print("TEMPLATE DNA LENGTH: ", len(template_dna))


		# Invertir o Conservar sentido de la secuencia para que attL quede a la izq #
		template_dna = checkSense(actual_orfx, template_dna)


		# Extraer 21 nucl finales de orfX que corresponden a attL #
		att_actual_orfx = actual_orfx[len(actual_orfx)-21:]
		print("ATT ACTUAL ORFX: ", att_actual_orfx)
		with open("attL_"+nombre+".fasta", "w") as f:
			f.write(">attL_"+nombre+"\n")
			f.write(att_actual_orfx+"\n")

		# Cortar template_dna desde orfx -> end
		texto_orfx, inicio_orfx, final_orfx = get_orfX_pos(template_dna, actual_orfx, "orfX")
		print("ORFX INFO: ", texto_orfx, inicio_orfx, final_orfx)
		template_dna = template_dna[inicio_orfx:]


		# ------------------------------------------------------------------------- #
		#            				Encontrar attR 				                    #

		args = get_attchment(nombre, raw_data, template_dna, attr_database_path, blastn_exe)

		if args is not None:
			print("\n"+"-"*78+"\n")
			print("ATTR ENCONTRADO")
			hit = args[2]

			texto_orfx, inicio_orfx, final_orfx = sequence_position(template_dna, actual_orfx, "orfX")
			print("\n"+"-"*78+"\n")
			print("COORDENADAS: ", inicio_orfx, args[1])
	
			sccmec = template_dna[inicio_orfx:args[1]]
			print("\n"+"-"*78+"\n")
			print "Contigs IDs: ", contig_ids
			print "Contig Length: ", len(template_dna)
			print "SCCmec Length: ", len(sccmec)
			print "attR match: ", hit
			print("\n"+"-"*78+"\n")


			texto_orfx, inicio_orfx, final_orfx = sequence_position(sccmec, actual_orfx, "orfX")
			print("ORFX INFO: ", texto_orfx, inicio_orfx, final_orfx)
			print("\n"+"-"*78+"\n")


			texto_mec, inicio_mec, final_mec = sequence_position(sccmec, actual_mecA, "PBP2a")
			print("PBP2a INFO: ", texto_mec, inicio_mec, final_mec)
			print("\n"+"-"*78+"\n")

			texto_ccr, inicio_ccr, final_ccr = sequence_position(sccmec, actual_ccr, "ccr")
			print("CCR INFO: ", texto_ccr, inicio_ccr, final_ccr)
			print("\n"+"-"*78+"\n")

			# ------------------------------------------------------------------------- #
			#        Crear archivo con la secuencia SCCmec desde attL hasta attR        #

			sccmec_file_name = "sccmec_"+nombre+".fasta" 
			with open(sccmec_file_name, "w") as f:
				f.write(">sccmec_"+nombre+"_l"+str(len(sccmec))+"\n")
				for i in range(0, len(sccmec), 60):
					f.write(sccmec[i:i+60]+'\n')

			# ------------------------------------------------------------------------- #
			#            Crear archivo con la informacion relevante                     #

			with open("info_"+nombre+".txt", "w") as f:
				f.write("Contig_ID:"+str(contig_ids)+'\n')
				f.write("Contig_length:"+str(len(template_dna))+'\n')
				f.write("SCCmec_length:"+str(len(sccmec))+'\n')
				f.write("attR_match:"+hit+'\n')
				f.write(texto_orfx+str(inicio_orfx)+'-'+str(final_orfx)+'\n')
				f.write(texto_mec+str(inicio_mec)+'-'+str(final_mec)+'\n')
				f.write(texto_ccr+str(inicio_ccr)+'-'+str(final_ccr)+'\n')
				f.write("Location:"+str(inicio_orfx)+"-"+str(args[1])+'\n')
				f.write("Quality:"+"A"+'\n')

			sys.exit()

		else:
			print("Deleting: ", raw_data)
			import shutil
			shutil.rmtree(raw_data)


	else:
		print("Different Contigs")
		os.chdir(raw_data)
		# ------------------------------------------------------------------------- #
		# ------------------------------------------------------------------------- #
		# ------------------------------------------------------------------------- #

		contig_left, contig_right = ordenar_contigs(contig_id_orfx, contig_id_mecA, contig_id_ccr)
		sequence_contig_left = get_sequence(contigs_dict, contig_left)
		sequence_contig_right = get_sequence(contigs_dict, contig_right)

		sequence_contig_left = checkSense(actual_orfx, sequence_contig_left)
		texto_orfx, inicio_orfx, final_orfx = sequence_position(sequence_contig_left, actual_orfx, "orfX")

		print("ORFX INFO: ", texto_orfx, inicio_orfx, final_orfx)

		
		sccmec_leftEnd = sequence_contig_left[inicio_orfx:]
		print("Largo SCCmec LeftEnd: ", len(sccmec_leftEnd))

		att_actual_orfx = actual_orfx[len(actual_orfx)-21:]
		print("ATT ACTUAL ORFX: ", att_actual_orfx)
		with open("attL_"+nombre+".fasta", "w") as f:
			f.write(">attL_"+nombre+"\n")
			f.write(att_actual_orfx+"\n")


		print("Largo Contig One: ", len(sequence_contig_left))
		print("Largo Contig Two: ", len(sequence_contig_right))

		# ------------------------------------------------------------------------- #
		# ------------------------------------------------------------------------- #
		# ------------------------------------------------------------------------- #

		args = get_attchment(nombre, raw_data, sequence_contig_right, attr_database_path, blastn_exe)

		if args is not None:
			print("\n"+"-"*78+"\n")
			print("ATTR ENCONTRADO")

			attR_inico, attR_final, hit = [x for x in args]

			sccmec_rightEnd = sequence_contig_right[:args[1]]

			sccmec = ''.join([sccmec_leftEnd, sccmec_rightEnd])

			print("\n"+"-"*78+"\n")
			print "Contigs IDs: ", contig_ids
			print "Contig Length: ", (len(sequence_contig_left)+len(sequence_contig_right))
			print "SCCmec Length: ", len(sccmec)
			print "attR match: ", hit
			print("\n"+"-"*78+"\n")


			texto_orfx, inicio_orfx, final_orfx = sequence_position(sccmec, actual_orfx, "orfX")
			print("ORFX INFO: ", texto_orfx, inicio_orfx, final_orfx)
			print("\n"+"-"*78+"\n")


			texto_mec, inicio_mec, final_mec = sequence_position(sccmec, actual_mecA, "PBP2a")
			print("PBP2a INFO: ", texto_mec, inicio_mec, final_mec)
			print("\n"+"-"*78+"\n")

			texto_ccr, inicio_ccr, final_ccr = sequence_position(sccmec, actual_ccr, "ccr")
			print("CCR INFO: ", texto_ccr, inicio_ccr, final_ccr)
			print("\n"+"-"*78+"\n")

			# ------------------------------------------------------------------------- #
			#        Crear archivo con la secuencia SCCmec desde attL hasta attR        #
			sccmec_file_name = "sccmec_"+nombre+".fasta" 
			with open(sccmec_file_name, "w") as f:
				f.write(">sccmec_"+nombre+"_l"+str(len(sccmec))+"\n")
				for i in range(0, len(sccmec), 60):
					f.write(sccmec[i:i+60]+'\n')

			# ------------------------------------------------------------------------- #
			#            Crear archivo con la informacion relevante                     #

			with open("info_"+nombre+".txt", "w") as f:
				f.write("Contig_ID:"+str(contig_ids)+'\n')
				f.write("Contig_length:"+str(len(sequence_contig_left)+len(sequence_contig_right))+'\n')
				f.write("SCCmec_length:"+str(len(sccmec))+'\n')
				f.write("attR_match:"+hit+'\n')
				f.write(texto_orfx+str(inicio_orfx)+'-'+str(final_orfx)+'\n')
				f.write(texto_mec+str(inicio_mec)+'-'+str(final_mec)+'\n')
				f.write(texto_ccr+str(inicio_ccr)+'-'+str(final_ccr)+'\n')
				f.write("Location:"+str(inicio_orfx)+"-"+str(len(sccmec))+'\n')
				f.write("Quality:"+"B"+'\n')

			sys.exit()

		sequence_contig_right = reverse_complement(sequence_contig_right)
		args = get_attchment(nombre, raw_data, sequence_contig_right, attr_database_path, blastn_exe)

		if args is not None:
			print("\n"+"-"*78+"\n")
			print("ATTR ENCONTRADO")

			attR_inico, attR_final, hit = [x for x in args]

			sccmec_rightEnd = sequence_contig_right[:args[1]]

			sccmec = ''.join([sccmec_leftEnd, sccmec_rightEnd])

			print("\n"+"-"*78+"\n")
			print "Contigs IDs: ", contig_ids
			print "Contig Length: ", (len(sequence_contig_left)+len(sequence_contig_right))
			print "SCCmec Length: ", len(sccmec)
			print "attR match: ", hit
			print("\n"+"-"*78+"\n")


			texto_orfx, inicio_orfx, final_orfx = sequence_position(sccmec, actual_orfx, "orfX")
			print("ORFX INFO: ", texto_orfx, inicio_orfx, final_orfx)
			print("\n"+"-"*78+"\n")


			texto_mec, inicio_mec, final_mec = sequence_position(sccmec, actual_mecA, "PBP2a")
			print("PBP2a INFO: ", texto_mec, inicio_mec, final_mec)
			print("\n"+"-"*78+"\n")

			texto_ccr, inicio_ccr, final_ccr = sequence_position(sccmec, actual_ccr, "ccr")
			print("CCR INFO: ", texto_ccr, inicio_ccr, final_ccr)
			print("\n"+"-"*78+"\n")

			# ------------------------------------------------------------------------- #
			#        Crear archivo con la secuencia SCCmec desde attL hasta attR        #

			sccmec_file_name = "sccmec_"+nombre+".fasta" 
			with open(sccmec_file_name, "w") as f:
				f.write(">sccmec_"+nombre+"_l"+str(len(sccmec))+"\n")
				for i in range(0, len(sccmec), 60):
					f.write(sccmec[i:i+60]+'\n')

			# ------------------------------------------------------------------------- #
			#            Crear archivo con la informacion relevante                     #

			with open("info_"+nombre+".txt", "w") as f:
				f.write("Contig_ID:"+str(contig_ids)+'\n')
				f.write("Contig_length:"+str(len(sequence_contig_left)+len(sequence_contig_right))+'\n')
				f.write("SCCmec_length:"+str(len(sccmec))+'\n')
				f.write("attR_match:"+hit+'\n')
				f.write(texto_orfx+str(inicio_orfx)+'-'+str(final_orfx)+'\n')
				f.write(texto_mec+str(inicio_mec)+'-'+str(final_mec)+'\n')
				f.write(texto_ccr+str(inicio_ccr)+'-'+str(final_ccr)+'\n')
				f.write("Location:"+str(inicio_orfx)+"-"+str(len(sccmec))+'\n')
				f.write("Quality:"+"B"+'\n')


		else:
			print("Picky Case")
			sys.exit()


		# ------------------------------------------------------------------------- #
		# ------------------------------------------------------------------------- #
		# ------------------------------------------------------------------------- #






if __name__ == '__main__':
	main()
