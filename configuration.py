from helpers import *

def main():
# ------------------------------------------------------------------------- #
#                   1. CONFIGURATION & OUTPUT SET UP                        #
# ------------------------------------------------------------------------- #
	import sys
	import os
	import ConfigParser

	config = ConfigParser.ConfigParser()

	config.readfp(open(r'config2.txt'))
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
	prefix = config.get('configuration', 'prefix_output')

	basedir = os.path.dirname(os.path.abspath(sys.argv[0]))

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


# ------------------------------------------------------------------------- #
#              2. LOOP THROUGH CONTIGS FILES IN INPUT                       #
# ------------------------------------------------------------------------- #

	contigs_to_analize = [f for f in os.listdir(input_contigs) if f.endswith('.fasta')]

	print os.getcwd()

	for contig in contigs_to_analize:
		print("\n"+"-"*78)*3
# ------------------------------------------------------------------------- #
#       a. output folder for contig prefix_<name>                           #

		contig_abspath = os.path.join(input_contigs, contig)
		filename = contig.split('.')[0]
		folder = prefix+filename
		output_folder = create_dir(global_output_dir, folder)
		print filename
		print output_folder

# ------------------------------------------------------------------------- #
#       IMPORTANTE
		os.chdir(output_folder)

# ------------------------------------------------------------------------- #
#       b. Run prokka & set up files                                        # 
		
		output_prokka = prefix+'prokka'
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
		pbp2a_hit, err = simpleBlast(blastp_exe, prot_db_path, pbp2a_base, 'PBP2a')
		pbp2a_hit = simpleBlastParser(pbp2a_hit)
		ccr_hit, err = simpleBlast(blastp_exe, prot_db_path, ccr_base, 'ccr')
		ccr_hit = simpleBlastParser(ccr_hit)

		print
		print('orfX Hit: ', orfx_nucl_hit)
		print('PBP2a Hit: ', pbp2a_hit)
		print('ccr Hit: ', ccr_hit)

		core_elements = [orfx_nucl_hit, pbp2a_hit, ccr_hit]

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
#       Check if CORE ELEMENTS are in the same contig                       #

		
		seq_region_id_orfx = get_region(gff, orfx_nucl_hit)
		seq_region_id_pbp2a = get_region(gff, pbp2a_hit)
		seq_region_id_ccr = get_region(gff, ccr_hit)
		region_ids = [seq_region_id_orfx, seq_region_id_pbp2a, seq_region_id_ccr]

		print('orfX in: ', seq_region_id_orfx)
		print('PBP2a in: ', seq_region_id_pbp2a)
		print('ccr in: ', seq_region_id_ccr)
		print('\n'+'-'*78+'\n')


		if check_region_id(region_ids):
			print region_ids





if __name__ == '__main__':
	main()