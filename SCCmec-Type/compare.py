import sys

type_I = ['orfX', 'IS431', 'mecA', 'mecR1-Truncated', 'IS1272', 'ccrB1-Truncated', 'ccrA1']
type_II = ['orfX', 'IS431', 'mecA', 'mecR1', 'MecI', 'ccrB2', 'ccrA2']
type_III = ['orfX', 'ccrC', 'IS431', 'IS431', 'IS431', 'IS431', 'mecA', 'mecR1-Truncated', 'MecI-Truncated', 'ccrB3', 'ccrA3']
type_IV = ['orfX', 'IS431', 'mecA', 'mecR1-Truncated', 'IS1272', 'ccrB2', 'ccrA2']
type_V = ['orfX', 'ccrC', 'IS431', 'mecA', 'ccrC']
type_VI = ['orfX-Truncated', 'IS431', 'mecA', 'mecR1-Truncated', 'IS1272', 'ccrB4', 'ccrA4']
type_VII = ['orfX', 'ccrC', 'IS431', 'mecA', 'mecR1-Truncated', 'IS431']
type_VIII = ['orfX', 'IS431', 'mecA', 'mecR1', 'MecI-Truncated', 'ccrB4', 'ccrA4']
type_IX = ['orfX', 'ccrB1', 'ccrA1', 'IS431-Truncated', 'mecA', 'IS431']
type_X = ['orfX', 'IS431', 'mecA', 'IS431', 'ccrB1-Truncated', 'ccrA1']
type_XI = ['orfX', 'mecC', 'mecR1', 'MecI', 'ccrB3', 'ccrA1']

tipos = [type_I, type_II, type_III, type_IV, type_V, type_VI, type_VII, type_VIII, type_IX, type_X, type_XI]

SCCmec_EU437549 = ['orfX', 'IS431', 'IS431', 'mecA', 'mecR1-Truncated', 'IS1272', 'ccrB2', 'ccrA2']


ccrA1B1 = ["ccrA1", "ccrB1"]
ccrA1B1_ = ["ccrA1", "ccrB1-Truncated"]
ccrA1B1__ = ["ccrA1-Truncated", "ccrB1"]

ccrA2B2 = ["ccrA2", "ccrB2"]
ccrA2B2_ = ["ccrA2", "ccrB2-Truncated"]
ccrA2B2__ = ["ccrA2-Truncated", "ccrB2"]

ccrA3B3 = ["ccrA3", "ccrB3"]
ccrA3B3_ = ["ccrA3", "ccrB3-Truncated"]
ccrA3B3__ = ["ccrA3-Truncated", "ccrB3"]

ccrA4B4 = ["ccrA4", "ccrB4"]
ccrA4B4_ = ["ccrA4", "ccrB4-Truncated"]
ccrA4B4__ = ["ccrA4-Truncated", "ccrB4"]

ccrC = ["ccrC"]

ccrA1B3 = ["ccrA1", "ccrB3"]

def ccrAllotype(lst):
	if all(gene in lst for gene in ccrA1B1):
		return "1"
	if all(gene in lst for gene in ccrA1B1_):
		return "1"	
	if all(gene in lst for gene in ccrA1B1__):
		return "1"

	if all(gene in lst for gene in ccrA2B2):
		return "2"
	if all(gene in lst for gene in ccrA2B2_):
		return "2"
	if all(gene in lst for gene in ccrA2B2__):
		return "2"

	if all(gene in lst for gene in ccrA3B3):
		return "3"
	if all(gene in lst for gene in ccrA3B3_):
		return "3"
	if all(gene in lst for gene in ccrA3B3__):
		return "3"

	if all(gene in lst for gene in ccrA4B4):
		return "4"
	if all(gene in lst for gene in ccrA4B4_):
		return "4"
	if all(gene in lst for gene in ccrA4B4__):
		return "4"

	if all(gene in lst for gene in ccrC):
		return "5"	

	if all(gene in lst for gene in ccrA1B3):
		return "8"

	else:
		return "N"


def mecClass2(lst):
	for i in range(len(lst)):
		if lst[i] == "mecA":

			if lst[i+1] == "mecR1":
				if lst[i+2] == "MecI":
					return "A"
			if lst[i+1] == "mecR1":
				if lst[i+2] == "MecI-Truncated":
					return "A"
			if lst[i+1] == "mecR1-Truncated":
				if lst[i+2] == "MecI":
					return "A"
			if lst[i+1] == "mecR1-Truncated":
				if lst[i+2] == "MecI-Truncated":
					return "A"

			if lst[i+1] == "mecR1-Truncated":
				if lst[i+2] == "IS1272":
					if lst[i-1] == "IS431":
						return "B"


			if lst[i-1] == "IS431":
				if lst[i+1] == "mecR1-Truncated":
					if lst[i+2] == "IS431":
						return "C1"

			if lst[i-1] == "IS431" or lst[i-1] == "IS431-Truncated":
				if lst[i+1] == "IS431" or lst[i+1] == "IS431-Truncated":

					try:
						if "ccrA" in lst[i-2] or "ccrB" in lst[i-2]:
							return "CL"
						if "ccrA" in lst[i+2] or "ccrB" in lst[i+2]:
							return "CR"
					except IndexError:
						pass

				return "C"


		if lst[i] == "mecC":
			if lst[i+1] == "mecR1":
				if lst[i+2] == "MecI":
					return "E"			




sccmec_dict = {
	"1B": "I",
	"2A": "II",
	"3A": "III",
	"2B": "IV",
	"5C": "V",
	"4B": "VI",
	"5C1": "VII",
	"4A": "VIII",
	"1CL": "IX",
	"1CR": "X",
	"8E": "XI"
}


for tipo in tipos:
	#print tipo
	clase = mecClass2(tipo)
	alotipo = ccrAllotype(tipo)
	code = (alotipo+clase)
	print sccmec_dict[code]

