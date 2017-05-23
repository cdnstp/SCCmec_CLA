import sys

data = []
db_dict_ident = {}
db_dict_align = {}

with open(sys.argv[1]) as f:
	lines = f.readlines()
	for line in lines:
		if not line:
			continue
		if line:
			args = line.split('\t')
			if args[0] != args[2]:
				pares = (sorted([args[0], args[2]]))
				pares = '-'.join([pares[0], pares[1]])
				#print pares+'\t'+args[10]

				porcentaje = (((float(args[5])-float(args[4]))+1)/float(args[1]))*100

				data.append(pares)

				if pares in db_dict_ident:
					db_dict_ident[pares].append(args[10])
				else:
					db_dict_ident[pares] = [args[10]]


				if pares in db_dict_align:
					db_dict_align[pares].append(porcentaje)
				else:
					db_dict_align[pares] = [porcentaje]


data = list(set(data))


for dat in data:
	query = db_dict_ident.pop(dat)
	if len(query) >= 2:
		#align = db_dict_align.pop(dat)
		query.sort(key=float)
		par1, par2 = dat.split('-')
		#print(par1+'\t'+par2+'\t'+query[0]+'\t'+str(align[0]))
		print(par1+" (ss) "+par2+" = "+query[0])
		#print(par1+" ss "+par2)


