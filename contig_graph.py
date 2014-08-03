
import networkx as nx

import matplotlib
matplotlib.use('Agg')

import Queue



def sum_of_all_nbrs_exept_longest(nbrs,contigs):
	lengths = sorted (map(lambda x: contigs[x[0]].length, nbrs))
	lengths.sort(reverse=True)
	return sum(lengths[1:])

def sum_of_repeat_lengths(path,contigs):
	assert len(path)%2 == 0
	structure_length = 0
	for repeat in path[::2]:
		structure_length += contigs[repeat[0]].length
	return structure_length

class ContigGraph(nx.MultiGraph):
	"""docstring for ContigGraph"""
	def __init__(self):
		super(ContigGraph, self).__init__()
		

	def __str__(self):
		graph_repr = ''
		# nodes_visited =set()
		# for edge in self.edges(data=True):
		# 	graph_repr += str(edge) +'\n'

		return(graph_repr)

	def copy_no_edge_info(self):
		G_copy = ContigGraph()
		#for edge in self.edges_iter():
			#print 'hgbfdv',edge
			#G_copy.add_edge(edge[0],edge[1])
		G_copy.add_edges_from([edge for edge in self.edges_iter() ])
		return G_copy


	def other_end(self,node):
		if node[1]:
			return (node[0],0)
		else:
			return (node[0],1)

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

	def add_link_edge(self,node1, orientation1, node2, orientation2, link_count, link_type, mean_obs,lib):
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
			Returns neighbors (from all libraries) of a node except itself as a neighbor.
			(If the contig node is linking to its other end)
		'''
		return filter(lambda x: x[0] != node[0],super(ContigGraph, self).neighbors(node))

	def lib_neighbors(self,node,lib):
		'''
			Returns neighbors (from a given library) of a node also including its
			other end as a neighbor. (If the contig node is linking to its other end)
		'''
		lib_nbrs = []
		for nbr in super(ContigGraph, self).neighbors(node):
			for lib_edge in self[node][nbr]:
				#print 'lib name:',self[node][nbr][lib_edge]['l'].lib_name, 'lib:', lib
				if self[node][nbr][lib_edge]['l'] == lib:
					lib_nbrs.append(nbr)
		return lib_nbrs

	def all_edges_between_two_nodes(self,node1,node2):
		edges =[]
		for lib_edge in self[node1][node2]:
			edges.append(self[node1][node2][lib_edge])
		return edges

	def effective_nr_of_links(self,lib_edge):
		if lib_edge['t'] == 'mp_am' or lib_edge['t'] == 'pe_am':
			return lib_edge['n'] * lib_edge['l'].contamine_rate
		else:
			return lib_edge['n']
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
						if sum_of_all_nbrs_exept_longest(self.lib_neighbors(node, lib), contigs) > lib.mean + 4*lib.sd + 2*besst.config_params.kmer_overlap and contigs[node[0]].coverage > 3 * besst.tot_mean_coverage:
							print 'here okay', sum_of_all_nbrs_exept_longest(self.lib_neighbors(node, lib), contigs)
							# for nbr in self.lib_neighbors(node, lib):
							# 	print 'info:',self[node][nbr]
							#print 'nbrs:',
							contigs[node[0]].is_repeat = True
							yield node
						else:
							pass
							#print 'THIS IS okay', sum_of_all_nbrs_exept_longest(self.lib_neighbors(node, lib), contigs)
					elif lib.lib_type == 'mp':
						if sum_of_all_nbrs_exept_longest(self.lib_neighbors(node, lib), contigs) > lib.mean + 4*lib.sd + 2*besst.config_params.kmer_overlap and contigs[node[0]].coverage > 3 * besst.tot_mean_coverage:
							print 'here okay', sum_of_all_nbrs_exept_longest(self.lib_neighbors(node, lib), contigs)
							contigs[node[0]].is_repeat = True
							yield node 
						else:
							pass
							#print 'THIS IS okay', sum_of_all_nbrs_exept_longest(self.lib_neighbors(node, lib), contigs)


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

	def get_contigs_flanking_a_repeat_structure(self, G_full, repeat_structure_left_end,repeat_structure_right_end):
		"""
			Parses the neighborhood region of one repeat R (or several repeats in a row e.g. R_1,R_2) 
			and looks for links spanning over it (them).
			If there are links between any two neighbors A and B from a given library, it will store
			the region ARB as a triplet. This triplet is found by serching for shortest paths between
			the two nodes representing R. THere might be several different triplets for a repeat. 
		"""
		try:
			for path in nx.all_shortest_paths(G_full, repeat_structure_left_end, repeat_structure_right_end):
				# return the flanking contig ends of a repeat (only if flanks are directly connected )
				if len(path) == 4:
					#print path[1],path[-2]
					# return flanting contig nodes plus the library object that is connecting them
					key = G_full[path[1]][path[-2]].keys()[0]
					yield path[1], path[-2], G_full[path[1]][path[-2]][key]['l']
			else:
				yield None,None, None
		except nx.exception.NetworkXNoPath:
			#print 'no path, i.e. repeat structure cannot be bridged over'
			yield None,None, None

	def repeat_structure_iterator(self, G_full, contigs,besst):
		"""
			Assemble complicated repeat regions using an iterative buildup of assembling the closest 
			neighbors of a repeat. This iterative procedure is continued until both endpoints are classified as
			non repeats. The "endpoints, as mentioned here is always connected in the original contig graph. Oherwise,
			we would have no information of joining the two endpoints.
		"""
		queue = Queue.Queue()
		for repeat in self.repeat_regions:
			queue.put([(repeat,0),(repeat,1)])


		while not queue.empty():
			repeat_structure = queue.get()

			repeat_structure_left_end,repeat_structure_right_end = repeat_structure[0], repeat_structure[-1]
			#for item in self.get_contigs_flanking_a_repeat_structure(G_full, repeat_structure_left_end, repeat_structure_right_end):
			#		print item

			for flank1,flank2,lib_object in self.get_contigs_flanking_a_repeat_structure(G_full, repeat_structure_left_end, repeat_structure_right_end):
				if flank1 == None and flank2 == None:
					# return unsuccessful repats here
					yield None
					continue
				if not contigs[flank1[0]].is_repeat and not contigs[flank2[0]].is_repeat:
					# bothe ends are not repeats
					# finally, check that the repeat structure is not spuriously build by
					# looking if at least one of the unique ends link to all of the repeat
					# contigs
					new_structure = [flank1] + repeat_structure + [flank2]
					yield  new_structure


				elif len(repeat_structure)/2 > 2:

					# heuristic stopping criterion
					continue

				elif contigs[flank1[0]].is_repeat and not contigs[flank2[0]].is_repeat:
					# left end is still a repeat, continue building structure
					#print '1'
					new_structure = [self.other_end(flank1), flank1] + repeat_structure 
					if sum_of_repeat_lengths(new_structure[:-2] ,contigs) <= lib_object.mean + 3*lib_object.sd + (len(new_structure[:-2])/2)*besst.config_params.kmer_overlap:
					
						queue.put(new_structure)

				elif not contigs[flank1[0]].is_repeat and contigs[flank2[0]].is_repeat:
					#right end is still a repeat, continue building structure
					#print '2'
					new_structure =  repeat_structure + [flank2,self.other_end(flank2)]
					if sum_of_repeat_lengths(new_structure[2:] ,contigs) <= lib_object.mean + 3*lib_object.sd + (len(new_structure[2:])/2)*besst.config_params.kmer_overlap:
						queue.put(new_structure)

				else:
					# Both ends are still repeats,  continue building structure
					#print '3'
					new_structure = [self.other_end(flank1), flank1] + repeat_structure + [flank2,self.other_end(flank2)]
					if sum_of_repeat_lengths(new_structure ,contigs) <= lib_object.mean + 3*lib_object.sd + (len(new_structure)/2)*besst.config_params.kmer_overlap:
						queue.put(new_structure)


	def most_supported_repeats(self,G_full, contigs,besst):
		repeat_structures = []
		for repeat_region in self.repeat_structure_iterator(G_full, contigs, besst):
			if not repeat_region:
				continue
			r = RepeatStructure(repeat_region)
			unique_flank1 = repeat_region[0]
			unique_flank2 = repeat_region[-1]

			if r.good_repeat_structure(G_full, unique_flank1, unique_flank2, repeat_region[1:-1]):
				#print 'lolsss'
				r.nr_supporting_links(G_full, unique_flank1, unique_flank2, repeat_region[1:-1])
				repeat_structures.append(r)
			else:
				pass#print 'spurious repeat structure'

		return sorted(repeat_structures, reverse=True)


class RepeatStructure(object):
	"""docstring for RepeatStructure"""
	def __init__(self, repeat_path):
		super(RepeatStructure, self).__init__()
		self.repeat_path = repeat_path

	def __lt__(self,other):
		if self.link_support < other.link_support:
			return True
		else:
			return False

	def __str__(self):
		return 'Repeat region with {0} supporting links. Number of links between unique contigs is {1}. Path is:\n{2}'.format(self.link_support,self.unique_ctgs_link_strength,self.repeat_path)
	def nr_supporting_links(self, G_full, unique_ctg1, unique_ctg2,repeat_structure):
		nr_left_links = map(lambda x: G_full[x][unique_ctg1][ G_full[x][unique_ctg1].keys()[0] ]['n'] if x in G_full.neighbors(unique_ctg1) else 0, repeat_structure[::2])
		nr_right_links = map(lambda x: G_full[x][unique_ctg2][ G_full[x][unique_ctg2].keys()[0] ]['n'] if x in G_full.neighbors(unique_ctg2) else 0, repeat_structure[1::2])
		tot_links = sum(nr_right_links)+sum(nr_left_links)
		self.link_support = tot_links
		self.unique_ctgs_link_strength = G_full[unique_ctg1][unique_ctg2][ G_full[unique_ctg1][unique_ctg2].keys()[0] ]['n']
		return tot_links	

	def good_repeat_structure(self, G_full, unique_ctg1, unique_ctg2,repeat_structure):
		left_links = map(lambda x: 1 if x in G_full.neighbors(unique_ctg1) else 0, repeat_structure[::2])
		right_links = map(lambda x: 1 if x in G_full.neighbors(unique_ctg2) else 0, repeat_structure[1::2])
		is_linked = []

		for l,r in zip(left_links, right_links):
			if l+r > 0:
				is_linked.append(True)
			else:
				is_linked.append(False)

		return(all(is_linked))	
