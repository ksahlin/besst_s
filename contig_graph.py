from networkx import MultiGraph


class ContigGraph(MultiGraph):
	"""docstring for ContigGraph"""
	def __init__(self):
		super(ContigGraph, self).__init__()
		

	def __str__(self):
		graph_repr = ''
		nodes_visited =set()
		for edge in self.edges(data=True):
			graph_repr += str(edge) +'\n'

		return(graph_repr)

	def add_node(self,node):
		""" True == Right end of sequence,
			False == Left end of sequence 
			Is sequence links to another sequence from 'True' end, 
			it is forward oriented, otherwise it is reverse complemented"""
		super(ContigGraph, self).add_nodes_from([(node,True),(node,False)])

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
			Returns neighbors of a node except itself as a neighbor. (If the contig node is linking to its other end)
		'''
		return filter(lambda x: x[0] != node[0],super(ContigGraph, self).neighbors(node))


	def remove_nbr_edges(self,node):
		for nbr in self.neighbors(node):
			self.remove_edge(node,nbr)

	def store_repeat_regions(self):
		"""
			remove a repeat and it's neighbors based on coverage and number of neighbors 
			We need to look at the number of unique neighboring contigs since its a multi graph. Thus,
			a there can be several edges between the two contigs
		"""
		self.repeat_regions = {}
		return
