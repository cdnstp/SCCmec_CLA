import sys

data = []
db_dict_ident = {}


with open(sys.argv[1]) as f:
	lines = f.readlines()
	for line in lines:
		if not line:
			continue
		if line:

			args = line.split('\t')

			sccmec_1 = args[0].split('/')[6].split('.')[0].replace('sccmec_', '')
			sccmec_2 = args[1].split('/')[6].split('.')[0].replace('sccmec_', '')


			if sccmec_1 != sccmec_2:

				pares = (sorted([sccmec_1, sccmec_2]))

				pares = '-'.join([pares[0], pares[1]])

				sim = "{:.6f}".format(1 - float(args[2]))

				#print sccmec_1, sccmec_2, sim

				data.append(pares)

				if pares in db_dict_ident:
					db_dict_ident[pares].append(sim)
				else:
					db_dict_ident[pares] = [sim]

#print len(data)
data = list(set(data))
#print len(data)


for dat in data:
	query = db_dict_ident.pop(dat)
	if len(query) >= 2:
		query.sort(key=float)
		par1, par2 = dat.split('-')

		#print(par1+" (ss) "+par2+" = "+query[0])
		print(par1+" ss "+par2)