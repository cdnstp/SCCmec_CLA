import sys
import networkx as nx


def input_network(file):
	""" toma archivo *.eda y lo transforma en una lista
		de listas con tres valores (nodo A, nodo B y similitud) """
	data = []
	with open(file) as f:
		lines = f.readlines()
		for line in lines:
			if not line:
				continue
			if line:
				source, _, target, _, weigth  = line.split()
				#source = source.split('_', 1)[1]
				#target = target.split('_', 1)[1]
				data.append([source, target, weigth])

	return data


def to_graph(rows):
	""" Agregar lista de lista a un objeto Graph de networkx
		nodo A, nodo B y similitud como weight de la conexion """
	G = nx.Graph()
	for row in rows:
		G.add_edge(row[0], row[1], weight=row[2])

	return G


def graph_cutoff(graph, threshold):
	""" Toma red inicial y la filtra segun valor de corte seleccionado,
	    devuelve un nuevo objeto Graph() con valores de edges >= T """
	filtered_graph = nx.Graph()
	for item in graph.edges(data=True):
		value = float(item[2]['weight'])
		if value >= threshold:
			filtered_graph.add_edge(str(item[0]), str(item[1]), weight = float(value))

	return filtered_graph		



def get_subgraph(graph, node):
	""" Entrega una lista de tuples con los pares de nodos
	del cluster en el que se encuentra un nodo dado """
	graphs = list(nx.connected_component_subgraphs(graph))

	for subgraph in graphs:
		for nodes in subgraph:
			if node in nodes:
				subgraph_con = subgraph.edges()

				return subgraph_con



def connected_components_network(edges, network):
	""" reconstruye una red de acuerdo a los pares de nodos dados """
	ccn = nx.Graph()

	for edge in edges:
		source, target = edge[0], edge[1]
		weight = network.get_edge_data(source, target, default = 0)
		value = weight['weight']

		ccn.add_edge(source, target, weight = float(value))

	return ccn


def cytoscape_files(graph, node):
	""" networkX Graph() -> Cytoscape datafiles """
	sif = '{0}_cytoscape_network.sif'.format(node)
	eda = '{0}_cytoscape_network.eda'.format(node)

	with open(eda, 'w') as f:
		f.write('Similarity (class=java.lang.Double)'+'\n')

	with open(sif, 'w') as f, open(eda, 'a') as g:
		for s,t,w in graph.edges(data=True):
			f.write(s+' ss '+t+'\n')
			g.write(s+' (ss) '+t+' = '+str(w['weight'])+'\n')




def related_cassettes(graph, node):
	""" ejecuta algoritmo de busqueda  """
	import operator
	sim, path = nx.single_source_dijkstra(graph, node)
	sorted_sim = sorted(sim.items(), key=operator.itemgetter(1))
	neighbors = graph.neighbors(node)

	file_name = '{0}_related_cassettes.txt'.format(node)
	neighbors_file = '{0}_neighbors_cassettes.txt'.format(node)

	with open(file_name, 'w') as f:
		f.write('ID\tDist(Sum of weight)\n')
		f.write('\n'.join('%s\t%s' % x for x in sorted_sim)  + '\n')

	with open(neighbors_file, 'w') as f:
		f.write(''.join('{0}\n'.format(x) for x in neighbors))		

	return sorted_sim




def draw_network(graph, node):
	""" opcional dado que calidad depende de la densidad """

	import matplotlib.pyplot as plt
	import matplotlib.patches as mpatches

	# HACK
	plt.close()
	# HACK
	
	# fijar color rojo para nodo (cassette) analizado y gris
	# para los demas

	color_map = []
	for nodes in graph:
		if nodes == node:
			color_map.append('red')
		else:
			color_map.append('grey')


	# color de las conexiones de acuerdo a similitud
	# normalizada al valor maximo
	weights = []
	for par in graph.edges():
		s, t = par[0], par[1]
		w = graph.get_edge_data(s, t)
		sim = w['weight']
		weights.append(sim)

	norm_weights = [float(n)/max(weights) for n in weights]

	# tamano de la figura
	plt.figure(1 , figsize = (20, 10))

	# layout de la red (depende de la densidad) 
	pos = nx.spring_layout(graph)

	nodes = nx.draw_networkx_nodes(graph, pos, 
								node_color = color_map,
								node_size = 1400,
								border = 12)
	edges = nx.draw_networkx_edges(graph, pos, 
								edge_color=norm_weights, 
								width=1, 
								edge_cmap=plt.cm.jet)
	labels = nx.draw_networkx_labels(graph, pos, font_size = 10)

	# colorbar para el valor de las conexiones
	cbar = plt.colorbar(edges)
	cbar.ax.set_ylabel('Normalized Edge Weight',
						labelpad=32, rotation=270, fontsize=12)
	
	# titulo
	fancy_title = 'Connected Components Graph of {0}'.format(node)
	plt.title(fancy_title, fontsize=16)

	# recuadro-margen
	frame = plt.gca()
	frame.axes.get_xaxis().set_ticks([])
	frame.axes.get_yaxis().set_ticks([])

	# leyenda, color de nodos
	node_circle = mpatches.Circle([0,0], radius = 1, color = 'red')
	related_circle = mpatches.Circle([0,0], radius = 1, color = 'grey')
	plt.legend([node_circle, related_circle], [node, 'related cassettes'], 
				loc='lower right', title="Node Color")

	# guardar figura
	plt.savefig('fig.png', dpi=300, orientation='landscape')
	plt.close()



def mash_exe(mash, chunk, sccmec_file):
	from subprocess import check_output
	output_mash = check_output([mash, 'dist', chunk, sccmec_file])
	return output_mash

#output_mash = mash_exe(mash_path, chunk, sys.argv[1])


def mash_parser(output_mash):
	data = []
	for line in output_mash.split('\n'):
		if not line:
			continue
		file_a, file_b, value, _, _ = line.split('\t')
		cassette_a = file_a.split('/')[-1].split('.')[0]
		cassette_b = file_b.split('/')[-1].split('.')[0]
		sim = '{:.6f}'.format(1 - float(value))
		data.append([cassette_a, cassette_b, sim])
	return data

#data = mash_parser(output_mash)


def network_classification(T, base_network, mash_path, chunk, sccmec_file, node):
# -------------------------------------------------------------------------- #
#   Setup graph database
	data = input_network(base_network)
	G = to_graph(data)
# -------------------------------------------------------------------------- #
#	Add new nodes/edges based on mash comparison 
	output_mash = mash_exe(mash_path, chunk, sccmec_file)
	data = mash_parser(output_mash)
	for row in data:
		G.add_edge(row[0], row[1], weight=row[2])

# -------------------------------------------------------------------------- #
#   Filter Graph based on current threshold or till find cluster
	threshold = float(T)
	subgraph_con = []
	while not subgraph_con:
		# test if subgraph_con exists at given threshold, 
		# otherwise select a less strict threshold.
		filtered_graph = graph_cutoff(G, threshold)
		subgraph_con = get_subgraph(filtered_graph, node)
		if not subgraph_con:
			threshold -= 0.001000

	else:
		if float(T) != threshold:
			exist = True
		else:
			exist = False
		ccn = connected_components_network(subgraph_con, G)
		cytoscape_files(ccn, node)
		draw_network(ccn, node)

		print len(ccn.nodes())

		sorted_sim = related_cassettes(ccn, node)

		return sorted_sim[1][0], exist
