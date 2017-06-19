def main():
	import sys
	import os
	import ConfigParser


	print sys.argv[0]
	base = os.path.dirname(os.path.abspath(sys.argv[0]))
	example = os.path.join(base, sys.argv[0])
	print example

	config = ConfigParser.ConfigParser()
	config.readfp(open(r'config2.txt'))
	# external programs
	prokka_exe = config.get('configuration', 'prokka')
	blastn_exe = config.get('configuration', 'blastn')
	blastp_exe = config.get('configuration', 'blastp')
	makeblastdb_exe = config.get('configuration', 'makeblastdb')
	# support data
	ccr = config.get('configuration', 'ccr')
	orfX = config.get('configuration', 'orfX')
	mecA = config.get('configuration', 'pbp2a')
	attr_db = config.get('configuration', 'attr_db')
	# inputs
	contigs_to_analize = config.get('configuration', 'input')
	# output dir
	output_data = config.get('configuration', 'output') 
	print mecA
	print base
	pbp2a_file = os.path.join(base, mecA)
	print pbp2a_file
	with open(pbp2a_file, 'r') as f:
		print f.readlines()
	pass


if __name__ == '__main__':
	main()