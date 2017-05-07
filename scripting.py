
def checkSense(sequence, contig):
	if re.findall("{pattern}".format(pattern=sequence), contig):
		return contig, "+"
	contig = reverse_complement(contig)
	if re.findall("{pattern}".format(pattern=sequence), contig):
		return contig, "+"


