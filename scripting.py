		if len(set(contigs_id)) == 2:
			if contig_id_mecA is contig_id_ccr:
				contigID_mec_ccr = contig_id_mecA
				print("mec and ccr in same contig: ", contigID_mec_ccr)
			if contig_id_mecA is contig_id_orfx:
				contigID_mec_orfx = contig_id_mecA
				print("mec and orfx in same contig: ", contigID_mec_orfx)
		else:
			print("I haven't reached that level of magnificence")
			sys.exit()

		template_dna_orfx = get_sequence(contigs_dict, contig_id_orfx)

		actual_orfx = get_sequence(nucl_dict, orfx_nucl_hit)

		#print actual_orfx
		texto_orfx, inicio_orfx, final_orfx = find_position(template_dna_orfx, actual_orfx, "orfX")

		print('\n'+"-"*78+'\n')
		print("Largo TEMPLATE DNA WITH ORFX: ", len(template_dna_orfx))
		print(texto_orfx, inicio_orfx, final_orfx)




		template_dna = get_sequence(contigs_dict, contig_id_mecA)
		actual_mecA = get_sequence(nucl_dict, mecA_hit)
		#print actual_mecA
		texto_mec, inicio_mec, final_mec = find_position(template_dna, actual_mecA, "PBP2a")
		print('\n'+"-"*78+'\n')
		print("Largo TEMPLATE DNA WITH MEC: ", len(template_dna))
		print(texto_mec, inicio_mec, final_mec)

		template_dna = get_sequence(contigs_dict, contig_id_ccr)
		actual_ccr = get_sequence(nucl_dict, ccr_hit)
		#print actual_mecA
		texto_ccr, inicio_ccr, final_ccr = find_position(template_dna, actual_ccr, "ccr")
		print('\n'+"-"*78+'\n')
		print("Largo TEMPLATE DNA WITH CCR: ", len(template_dna))
		print(texto_ccr, inicio_ccr, final_ccr)


		print('\n'+"-"*78+'\n')




def find_position(sequence, query, name):
	if re.findall("{pattern}".format(pattern=query), sequence):
		for i in re.findall("{pattern}".format(pattern=query), sequence):
			print "+_%s_located_at:" % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]
			x = "+_%s_located_at:" % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]
			return x
	else:
		for i in re.findall("{pattern}".format(pattern=reverse_complement(query)), sequence):
			print "-_%s_located_at:" % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]
			x = "-_%s_located_at:" % name, [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), sequence)]
			return x




		template_dna_orfx = get_sequence(contigs_dict, contig_id_orfx)

		print actual_orfx
		texto_orfx, inicio_orfx, final_orfx = find_position(template_dna_orfx, actual_orfx, "orfX")
		print('\n'+"-"*78+'\n')
		print("Largo TEMPLATE DNA WITH ORFX: ", len(template_dna_orfx))
		print(texto_orfx, inicio_orfx, final_orfx)

		template_dna = get_sequence(contigs_dict, contig_id_mecA)
		actual_mecA = get_sequence(nucl_dict, mecA_hit)
		#print actual_mecA
		texto_mec, inicio_mec, final_mec = find_position(template_dna, actual_mecA, "PBP2a")
		print('\n'+"-"*78+'\n')
		print("Largo TEMPLATE DNA WITH MEC: ", len(template_dna))
		print(texto_mec, inicio_mec, final_mec)

		template_dna = get_sequence(contigs_dict, contig_id_ccr)
		actual_ccr = get_sequence(nucl_dict, ccr_hit)
		#print actual_mecA
		texto_ccr, inicio_ccr, final_ccr = find_position(template_dna, actual_ccr, "ccr")
		print('\n'+"-"*78+'\n')
		print("Largo TEMPLATE DNA WITH CCR: ", len(template_dna))
		print(texto_ccr, inicio_ccr, final_ccr)

		print('\n'+"-"*78+'\n')

		