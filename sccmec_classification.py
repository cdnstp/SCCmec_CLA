from helpers import *

def main():
# ------------------------------------------------------------------------- #
#                   1. CONFIGURATION & OUTPUT SET UP                        #
# ------------------------------------------------------------------------- #
	import sys
	import os
	import ConfigParser

	config = ConfigParser.ConfigParser()

	config.readfp(open(r'config.txt'))
	prokka = config.get('configuration', 'prokka')
	blastn = config.get('configuration', 'blastn')
	blastp = config.get('configuration', 'blastp')
	makedb = config.get('configuration', 'makeblastdb')
	ccr = config.get('configuration', 'ccr')
	orfX = config.get('configuration', 'orfX')
	pbp2a = config.get('configuration', 'pbp2a')
	attr_db = config.get('configuration', 'attr_db')
	input_dir = config.get('configuration', 'input_dir')
	output_data = config.get('configuration', 'output')
	threshold = config.get('configuration', 'threshold')
	core_db = config.get('configuration', 'core_db')
	base_network = config.get('configuration', 'base_network')
	chunk = config.get('configuration', 'chunk')
	mash = config.get('configuration', 'mash')

	basedir = os.path.dirname(os.path.abspath(sys.argv[0]))

	mash_path = os.path.join(basedir, mash)
	makeblastdb_exe = os.path.join(basedir, makedb)
	prokka_exe = blastp_exe = os.path.join(basedir, prokka)
	blastn_exe = os.path.join(basedir, blastn)
	blastp_exe = os.path.join(basedir, blastp)
	input_contigs = os.path.join(basedir, input_dir)
	orfx_base = os.path.join(basedir, orfX)
	ccr_base = os.path.join(basedir, ccr)
	pbp2a_base = os.path.join(basedir, pbp2a)
	global_output_dir = create_dir(basedir, output_data)
	attr_database_path = os.path.join(basedir, attr_db)
	core_db_path = os.path.join(basedir, core_db)
	base_network_path = os.path.join(basedir, base_network)
	chunk_path = os.path.join(basedir, chunk)



# ------------------------------------------------------------------------- #
#              2. LOOP THROUGH CONTIGS FILES IN INPUT                       #
# ------------------------------------------------------------------------- #

	contigs_to_analize = [f for f in os.listdir(input_contigs) if f.endswith('.fasta')]

	print os.getcwd()

	for contig in contigs_to_analize:
		print("\n"+"-"*78)*3
# ------------------------------------------------------------------------- #
#       a. output folder for contig output_<name>                           #

		contig_abspath = os.path.join(input_contigs, contig)
		filename = contig.split('.')[0]
		folder = 'output_{0}'.format(filename)
		output_folder = create_dir(global_output_dir, folder)
		print filename
		print output_folder

# ------------------------------------------------------------------------- #
		os.chdir(output_folder)
# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
#       b. Run prokka & set up files                                        # 
		
		output_prokka = 'output_prokka_{0}'.format(filename)
		execute_prokka(prokka_exe, output_prokka, contig_abspath, filename)
		ffn, gff, fna, faa = prokka_files(output_prokka)
		ffn, gff, fna, faa = list(os.path.join(output_folder, x) for x in [ffn, gff, fna, faa])
		fna_dict = fasta2dict(fna)
		ffn_dict = fasta2dict(ffn)

# ------------------------------------------------------------------------- #
#       c. BLAST Databases configuration                                    # 

		nucl_db_dir = create_dir(output_folder, 'nucl_db_dir')
		nucl_db = makeblastdb(makeblastdb_exe, nucl_db_dir, ffn, 'nucl', 'nucl_db')
		nucl_db_path = os.path.join(nucl_db_dir, nucl_db)

		prot_db_dir = create_dir(output_folder, 'prot_db_dir')
		prot_db = makeblastdb(makeblastdb_exe, prot_db_dir, faa, 'prot', 'prot_db')
		prot_db_path = os.path.join(prot_db_dir, prot_db)
		print
# ------------------------------------------------------------------------- #
#       d. Check if is MRSA by running BLAST searches of core elements      #

		orfx_nucl_hit, err = simpleBlast(blastn_exe, nucl_db_path, orfx_base, 'orfX')
		orfx_nucl_hit = simpleBlastParser(orfx_nucl_hit)
		mec_hit, err = simpleBlast(blastp_exe, prot_db_path, pbp2a_base, 'PBP2a')
		mec_hit = simpleBlastParser(mec_hit)
		ccr_hit, err = simpleBlast(blastp_exe, prot_db_path, ccr_base, 'ccr')
		ccr_hit = simpleBlastParser(ccr_hit)

		print
		print('orfX Hit: ', orfx_nucl_hit)
		print('PBP2a Hit: ', mec_hit)
		print('ccr Hit: ', ccr_hit)

		core_elements = [orfx_nucl_hit, mec_hit, ccr_hit]

		if not all(core_elements):
			print('It\'s not MRSA')
			os.chdir(output_folder)
			info_file = 'MSSA_{0}.txt'.format(filename)
			with open(info_file, "w") as f:
				f.write('Gene\tLocus\n')
				f.write('orfX\t{0}\n'.format(str(core_elements[0])))
				f.write('PBP2a\t{0}\n'.format(str(core_elements[1])))
				f.write('ccr\t{0}\n'.format(str(core_elements[2])))
			if core_elements[0] is not None:
				actual_orfx = get_sequence(ffn_dict, orfx_nucl_hit)
				attB_sequence = actual_orfx[len(actual_orfx)-60:]
				attB_datafile = 'attB_{0}.fasta'.format(filename)
				with open(attB_datafile, "w") as f:
					f.write('>attB_{0}\n'.format(filename))
					f.write(attB_sequence+'\n')

			continue

		print
		print('It\'s MRSA'+'\n')
# ------------------------------------------------------------------------- #
		os.chdir(output_folder)
# ------------------------------------------------------------------------- #

		actual_orfx = get_sequence(ffn_dict, orfx_nucl_hit)
		actual_mec = get_sequence(ffn_dict, mec_hit)
		actual_ccr = get_sequence(ffn_dict, ccr_hit)



# ------------------------------------------------------------------------- #
#       e. Check if CORE ELEMENTS are in the same contig                    #

		seq_region_id_orfx = get_region(gff, orfx_nucl_hit)
		seq_region_id_pbp2a = get_region(gff, mec_hit)
		seq_region_id_ccr = get_region(gff, ccr_hit)

		region_ids = [seq_region_id_orfx, seq_region_id_pbp2a, seq_region_id_ccr]

		print('orfX in: ', seq_region_id_orfx)
		print('PBP2a in: ', seq_region_id_pbp2a)
		print('ccr in: ', seq_region_id_ccr)
		print('\n'+'-'*78+'\n')

		if check_region_id(region_ids):
# ------------------------------------------------------------------------- #
#       f. If all core elements in the same fna sequence                    #

			# Invertir o Conservar sentido de la secuencia para que attL quede a la izq #
			fna_sequence = get_sequence(fna_dict, seq_region_id_orfx)
			fna_sequence = checkSense(actual_orfx, fna_sequence)
			print('Contig Length', len(fna_sequence))


			# Continuar utilizando fna_sequence desde orfx
			texto_orfx, inicio_orfx, final_orfx = get_orfX_pos(fna_sequence, actual_orfx, "orfX")
			print("ORFX INFO: ", texto_orfx, inicio_orfx, final_orfx)

			template_dna = fna_sequence[inicio_orfx:]

# ------------------------------------------------------------------------- #
#       g. Encontrar attR                                                   #
			params = get_attchment(filename, 
								output_folder, 
								template_dna, 
								attr_database_path, 
								blastn_exe)

			if params:
				attr_s, attr_e, hit = params
				orfx_text, orfx_s, orfx_e = sequence_position(template_dna, actual_orfx, "orfX")
				# Extraer y Guardar attL
				att_actual_orfx = actual_orfx[len(actual_orfx)-60:]
				attL_datafile = 'attL_{0}.fasta'.format(filename)
				with open(attL_datafile, 'w') as f:
					f.write('>attL_{0}_{1}_{2}\n'.format(filename, orfx_s, orfx_e))
					f.write(att_actual_orfx+'\n')

				# cassette from orfX start to attR end
				cassette = template_dna[orfx_s:attr_e]

				### MUST FIX THIS TO AVOID FUTURE BUGS ###
				mec_text, mec_s, mec_e = sequence_position(cassette, actual_mec, "PBP2a")
				ccr_text, ccr_s, ccr_e = sequence_position(cassette, actual_ccr, "ccr")


# ------------------------------------------------------------------------- #
#       h. Crear archivo con la secuencia SCCmec desde attL hasta attR      #
#		   Listo para anotacion y clasificacion                             #
# ------------------------------------------------------------------------- #

				cassette_filename = 'sccmec_{0}.fasta'.format(filename) 
				with open(cassette_filename, 'w') as f:
					f.write('>sccmec_{0}_l{1}\n'.format(filename, str(len(cassette))))
					for i in range(0, len(cassette), 60):
						f.write(cassette[i:i+60]+'\n')

				info_file = 'MRSA_{0}.txt'.format(filename)
				with open(info_file, 'w') as f:
					f.write('fna length\t{0}\n'.format(str(len(template_dna))))
					f.write('SCCmec length\t{0}\n'.format(str(len(cassette))))				
					f.write('Quality\tA\n')

				cassette_path = os.path.join(output_folder, cassette_filename)

				
# ------------------------------------------------------------------------- #
#    	i. Si no encuentra un attR en el mismo fna en el que estan los core #
#          elements se detiene y pasa al siguiente contig en INPUT          #
# ------------------------------------------------------------------------- #
			else:
				info_file = 'attr_not_found_{0}.txt'.format(filename)
				with open(info_file, 'w') as f:
					f.write('attr not found\n')
				continue
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #

		else:
			print('core elements contained in more than one fna')
			fnas = sort_fna(seq_region_id_orfx, seq_region_id_pbp2a, seq_region_id_ccr)

			if len(fnas) == 3: 
				info_file = 'separated_core_elements_{0}.txt'.format(filename)
				with open(info_file, 'w') as f:
					f.write('each core elements in differents fna\n')
				continue

			else:
				os.chdir(output_folder)

				contig_left, contig_right = fnas
				sequence_contig_left = get_sequence(fna_dict, contig_left)
				sequence_contig_right = get_sequence(fna_dict, contig_right)

				sequence_contig_left = checkSense(actual_orfx, sequence_contig_left)
				orfx_text, orfx_s, orfx_e = sequence_position(sequence_contig_left, actual_orfx, "orfX")
				
				# Extraer y Guardar attL
				att_actual_orfx = actual_orfx[len(actual_orfx)-60:]
				attL_datafile = 'attL_{0}.fasta'.format(filename)
				with open(attL_datafile, 'w') as f:
					f.write('>attL_{0}_{1}_{2}\n'.format(filename, orfx_s, orfx_e))
					f.write(att_actual_orfx+'\n')

				cassette_left_end = sequence_contig_left[orfx_s:]

				params = get_attchment(filename, 
									output_folder, 
									sequence_contig_right, 
									attr_database_path, 
									blastn_exe)
				if params:
					attr_s, attr_e, hit = params
					cassette_right_end = sequence_contig_right[:attr_e]

					cassette = ''.join([cassette_left_end, cassette_right_end])
					### MUST FIX THIS TO AVOID FUTURE BUGS ###
					mec_text, mec_s, mec_e = sequence_position(cassette, actual_mec, "PBP2a")
					ccr_text, ccr_s, ccr_e = sequence_position(cassette, actual_ccr, "ccr")


					cassette_filename = 'sccmec_{0}.fasta'.format(filename) 
					with open(cassette_filename, 'w') as f:
						f.write('>sccmec_{0}_l{1}\n'.format(filename, str(len(cassette))))
						for i in range(0, len(cassette), 60):
							f.write(cassette[i:i+60]+'\n')

					info_file = 'MRSA_{0}.txt'.format(filename)
					with open(info_file, 'w') as f:
						f.write('fna-L length\t{0}\n'.format(str(len(sequence_contig_left))))
						f.write('fna-R length\t{0}\n'.format(str(len(sequence_contig_right))))
						f.write('SCCmec length\t{0}\n'.format(str(len(cassette))))				
						f.write('Quality\tB\n')

					cassette_path = os.path.join(output_folder, cassette_filename)

				else:
					sequence_contig_right = reverse_complement(sequence_contig_right)
					params = get_attchment(filename, 
										output_folder, 
										sequence_contig_right, 
										attr_database_path, 
										blastn_exe)
					if params:
						attr_s, attr_e, hit = params
						cassette_right_end = sequence_contig_right[:attr_e]

						cassette = ''.join([cassette_left_end, cassette_right_end])
						### MUST FIX THIS TO AVOID FUTURE BUGS ###
						mec_text, mec_s, mec_e = sequence_position(cassette, actual_mec, "PBP2a")
						ccr_text, ccr_s, ccr_e = sequence_position(cassette, actual_ccr, "ccr")


						cassette_filename = 'sccmec_{0}.fasta'.format(filename) 
						with open(cassette_filename, 'w') as f:
							f.write('>sccmec_{0}_l{1}\n'.format(filename, str(len(cassette))))
							for i in range(0, len(cassette), 60):
								f.write(cassette[i:i+60]+'\n')

						info_file = 'MRSA_{0}.txt'.format(filename)
						with open(info_file, 'w') as f:
							f.write('fna-L length\t{0}\n'.format(str(len(sequence_contig_left))))
							f.write('fna-R length\t{0}\n'.format(str(len(sequence_contig_right))))
							f.write('SCCmec length\t{0}\n'.format(str(len(cassette))))				
							f.write('Quality\tC\n')

						cassette_path = os.path.join(output_folder, cassette_filename)
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
#       i'. Si no encuentra attR al buscar en ambos sentidos (5'-3' o 3'-5')#
#           se detiene y pasa al siguiente contig en INPUT                  #
					else:
						info_file = 'attr_not_found_{0}.txt'.format(filename)
						with open(info_file, 'w') as f:
							f.write('attr not found\n')
						continue
# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
# ------------------------------------------------------------------------- #
#       Luego de extraer el cassette exitosamente continua con la           #
#       clasificacion segun estandar (IWG-SCC) y nuestro clasificacion      #


		print cassette_path

		print('\n{0} annotation using current methods\n'.format(filename))

		sccmec_id = current_annotation(cassette_path, blastp_exe, prokka_exe, core_db_path, output_folder)

		print base_network_path

		print('\n{0} annotation using network classification\n'.format(filename))
		
		network_classification(threshold, base_network_path, mash_path, chunk_path, cassette_path, sccmec_id)


	sys.exit('done')

if __name__ == '__main__':
	main()

