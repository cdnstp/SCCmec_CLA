def find_position(sequence, gene, name):
	if re.findall("{pattern}".format(pattern=gene), sequence):
		for i in re.findall("{pattern}".format(pattern=gene), sequence):
			print "+", "%s located at: ", [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), seq)] % (name)
	else:
		for i in re.findall("{pattern}".format(pattern=reverse_complementary(gene)), sequence):
			print "-", "%s located at: ", [(m.start(0), m.end(0)) for m in re.finditer("{}".format(i), seq)] % (name)
