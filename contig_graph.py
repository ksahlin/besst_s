import os

import networkx as nx

import matplotlib
matplotlib.use('Agg')


def sum_of_all_nbrs_exept_longest(nbrs,contigs):
	lengths = sorted (map(lambda x: contigs[x[0]].length, nbrs))
	lengths.sort(reverse=True)
	return sum(lengths[1:])


class ContigGraph(nx.MultiGraph):
	"""docstring for ContigGraph"""
	def __init__(self):
		super(ContigGraph, self).__init__()
		

	def __str__(self):
		graph_repr = ''
		nodes_visited =set()
		for edge in self.edges(data=True):
			graph_repr += str(edge) +'\n'

		return(graph_repr)

	def other_end(self,node):
		if node[1] == 0:
			return (node[0],1)
		else:
			return (node[0],0)

	def draw(self, path, nodes):
		nx.draw(self,nodelist=nodes)
		matplotlib.pyplot.savefig(path)
		matplotlib.pyplot.clf()
		#nx.write_dot(self,path) #super(ContigGraph, self).draw()

	def add_node(self,node):
		""" True == Right end of sequence,
			False == Left end of sequence 
			Is sequence links to another sequence from 'True' end, 
			it is forward oriented, otherwise it is reverse complemented"""
		super(ContigGraph, self).add_nodes_from([(node,True),(node,False)])

	def remove_node(self,node):
		""" True == Right end of sequence,
			False == Left end of sequence 
			Is sequence links to another sequence from 'True' end, 
			it is forward oriented, otherwise it is reverse complemented"""
		super(ContigGraph, self).remove_nodes_from([(node,True),(node,False)])	

	def remove_nodes_from(self,nodes):
		""" True == Right end of sequence,
			False == Left end of sequence 
			Is sequence links to another sequence from 'True' end, 
			it is forward oriented, otherwise it is reverse complemented"""
		for node in nodes:
			super(ContigGraph, self).remove_nodes_from([(node,True),(node,False)])

	def add_edge(self,node1, orientation1, node2, orientation2, link_count, link_type, mean_obs,lib):
		# #print self.neighbors((node2,orientation2))
		# #print self[node1] #dir(super(ContigGraph, self))
		# if  (node2,orientation2) in self and (node1,orientation1) in self.neighbors((node2,orientation2)):
		# 	print self[(node2,orientation2)][(node1,orientation1)]
		#else:
		super(ContigGraph, self).add_edge((node1,orientation1),(node2,orientation2), n = link_count, t=link_type, o_bar = mean_obs, l=lib ) 


	def remove_self_links(self):
		for edge in super(ContigGraph, self).edges():
			if edge[0][0] == edge[1][0]:
				self.remove_edge(edge[0],edge[1])


	def iteredges(self):
		'''
			Iterator over edges that links two different sequences (no self-edges connecting two ends of the same sequence)
		'''
		for edge in super(ContigGraph, self).edges_iter():
			if edge[0][0] != edge[1][0]:
				yield edge

	def neighbors(self,node):
		'''
			Returns neighbors (from all libraries) of a node except itself as a neighbor. (If the contig node is linking to its other end)
		'''
		return filter(lambda x: x[0] != node[0],super(ContigGraph, self).neighbors(node))

	def lib_neighbors(self,node,lib):
		'''
			Returns neighbors (from all libraries) of a node except itself as a neighbor. (If the contig node is linking to its other end)
		'''
		lib_nbrs = []
		for nbr in super(ContigGraph, self).neighbors(node):
			for lib_edge in self[node][nbr]:
				#print 'lib name:',self[node][nbr][lib_edge]['l'].lib_name, 'lib:', lib
				if self[node][nbr][lib_edge]['l'] == lib:
					lib_nbrs.append(nbr)
		return lib_nbrs
			#print self[node][nbr]
		#return filter(lambda x: x[0] != node[0] and self[node][x]['l'].lib_name == lib, super(ContigGraph, self).neighbors(node))

	# def remove_nbr_edges(self,node):
	# 	for nbr in self.neighbors(node):
	# 		self.remove_edge(node,nbr)

	def repeat_iter(self,contigs,besst):
		"""
			If a node has more than two neighbors on a side and the two longest of these neighbors
			are longer than lib_mean + 4*lib_std + 2*kmer_overlap, then we know that at least one link 
			must be wrong or it is a repeat (or a split haplotype region).
		"""
		for node in self.nodes():
			for lib in besst.libs: 
				if len(self.lib_neighbors(node, lib)) > 1:
					if lib.lib_type == 'pe':
						if sum_of_all_nbrs_exept_longest(self.lib_neighbors(node, lib), contigs) > lib.mean_innies + 4*lib.sd_innies + 2*besst.config_params.kmer_overlap:
							
							#Only small repeats in this data set? look for coverage?
							yield node
					elif lib.lib_type == 'mp':
						if sum_of_all_nbrs_exept_longest(self.lib_neighbors(node, lib), contigs) > lib.mean_outies + 4*lib.sd_outies + 2*besst.config_params.kmer_overlap:
							yield node 

	def freeze_repeat_regions(self, contigs, besst):
		"""
			remove a repeat and it's neighbors based on coverage and number of neighbors 
			We need to look at the number of unique neighboring contigs since its a multi graph. Thus,
			a there can be several edges between the two contigs
		"""
		self.repeat_regions = {}
		for repeat in self.repeat_iter(contigs,besst):
			if repeat[0] in self.repeat_regions:
				# We have already added the repeat because vi have visitid the
				# other end of it
				pass
			else:
				self.repeat_regions[repeat[0]] = self.subgraph([repeat,self.other_end(repeat)]+ self.neighbors(repeat) + self.neighbors(self.other_end(repeat)) )
			
			#repeats.append(repeat)

			print repeat[0], contigs[repeat[0]].coverage #G[repeat]
			#graph_path = os.path.join(besst.config_params.output_path,'region-{0}.png'.format(repeat[0]))	
			#self.repeat_regions[repeat[0]].draw(graph_path, self.repeat_regions[repeat[0]].nodes())

		self.remove_nodes_from([repeat for repeat in self.repeat_regions])

