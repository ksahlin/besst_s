import argparse
import sys
import os
import pickle

from random import choice

from parser import fasta_parser,config_parser,link_parser
import sequence
import read_library
from contig_graph import ContigGraph 
from assembly_metrics import AssemblyMetrics
import paths

class GlobaInfo(object):
	"""docstring for GlobaInfo"""
	def __init__(self):
		super(GlobaInfo, self).__init__()
		self.config_params = None
		self.libs = []

	def __str__(self):
		return

	def get_cut_vertex_cutoff(self):
		cutoff = 0
		for lib in self.libs:
			if lib.lib_type == 'mp':
				if cutoff <= lib.mean_outies + 3*lib.sd_outies:
					cutoff = lib.mean_outies + 3*lib.sd_outies

			elif lib.lib_type == 'pe':
				if cutoff <= lib.mean_innies + 3*lib.sd_innies:
					cutoff = lib.mean_innies + 3*lib.sd_innies	
		self.cutoff = cutoff
	def calculate_total_mean_coverage(self):
		tot_cov = 0
		for lib in self.libs:
			tot_cov += lib.mean_coverage
		self.tot_mean_coverage = tot_cov		

besst = GlobaInfo()

def read_in_contigs(config_params,contig_dict):
	for accession,seq in fasta_parser.get_contigs(config_params.contig_file):
		c = sequence.Contig(accession, sequence = seq)
		contig_dict[c.name] = c

def read_library_metrics(config_params, contigs):
	lib_name =1
	for  lib_info in config_params.libs:
		if len(lib_info) == 5:
			lib_type, aligner, lib_loc ,mean, sd = lib_info
			lib = read_library.Library(lib_name, lib_type, aligner, lib_loc, os.path.join(config_params.output_path, 'lib_{0}_links.txt'.format(lib_name)),mean=mean,sd=sd)

		elif len(lib_info) == 3:
			lib_type, aligner, lib_loc = lib_info
			lib = read_library.Library(lib_name, lib_type, aligner, lib_loc, os.path.join(config_params.output_path, 'lib_{0}_links.txt'.format(lib_name)))


		elif len(lib_info) == 6:
			lib_type, aligner, lib_loc1,lib_loc2 ,mean, sd = lib_info
			lib = read_library.Library(lib_name, lib_type, aligner, lib_loc1, os.path.join(config_params.output_path, 'lib_{0}_links.txt'.format(lib_name)), bam_path2=lib_loc2,mean=mean,sd=sd)
		
		elif len(lib_info) == 4:
			lib_type, aligner, lib_loc1,lib_loc2 = lib_info
			lib = read_library.Library(lib_name, lib_type, aligner, lib_loc1, os.path.join(config_params.output_path, 'lib_{0}_links.txt'.format(lib_name)), bam_path2=lib_loc2)


		lib.get_libary_metrics()
		#print contigs
		lib.calculate_coverage(contigs)
		lib_name += 1
		besst.libs.append(lib)
		print lib

		

def create_link_file(contigs):
	for lib in besst.libs:
		lib.sorted_bam_to_link_file(contigs,besst.config_params.kmer_overlap)


def create_graph(G,contigs):

	##
	# create nodes

	for ctg in contigs:
		G.add_node(contigs[ctg].name)
	##
	# Create edges
	for lib in besst.libs:
		for (ctg1, orientation1, ctg2, orientation2, link_count, link_type, mean_obs) in  link_parser.get_links(lib.link_file_path):
			if link_count >= besst.config_params.min_links:
				G.add_link_edge(ctg1, orientation1, ctg2, orientation2, link_count, link_type, mean_obs, lib)

				#G.add_edge((subsequences[ seq1 ] ,orientation1),(subsequences[ seq2 ],orientation2),d=naive_gap,s=score.nr_links(link_count))
	#print G
	G.remove_self_links()

# 	return(G)

def score_scaffolding(path_factory, G_full):
	good_links = 0
	bad_links = 0
	for path in path_factory.path_ends.itervalues():
		path.find_supporting_links(G_full)
		good_links += path.good_links
		bad_links += path.bad_links

	return good_links/float(bad_links)


def main(args):
	params = config_parser.ConfigParams()
	params.read_cfg(args.config_file)
	besst.config_params = params
	contigs = {}
	read_in_contigs(params,contigs)
	metrics = AssemblyMetrics(contigs)
	print metrics
	read_library_metrics(params,contigs)

	# for ctg in contigs:
	# 	print contigs[ctg].name
	# 	print contigs[ctg].length, contigs[ctg].coverage

	create_link_file(contigs)
	G = ContigGraph()
	create_graph(G,contigs)
	besst.calculate_total_mean_coverage()


	#pickle.dump( G.copy_no_edge_info(), open( "/tmp/save_graph.p", "wb" ) )
	pickle.dump( G, open( "/tmp/save_graph.p", "wb" ) )

	# filter out repeats temporarily

	G.freeze_repeat_regions(contigs, besst)


	# graph_path = os.path.join(besst.config_params.output_path,'ctg_graph_w_repeats.png')	
	# G.draw(graph_path,repeats)
	graph_path = os.path.join(besst.config_params.output_path,'ctg_graph_no_repeats1.png')
	# G.remove_nodes_from([repeat[0] for repeat in repeats])
	
	#G.draw(graph_path,G.nodes())
	print len(G.edges()),G.size(),len(G.nodes())

	#print G[('50__len__471', 1)][ ('82__len__373', 0)]
	#print G.edges(data=True)
	#print G.all_edges_between_two_nodes(('27__len__157', False), ('18__len__153', 1))




	# output paths in sorted score order
	G_full = pickle.load( open( "/tmp/save_graph.p", "rb" ) )
	#print G_full.edges(data=True)

	repeat_paths = {}
	tmp_file = open('/tmp/repeats.txt', 'w')
	unique_flanks = set()

	for repeat_structure in G.most_supported_repeats(G_full, contigs,besst):
		if repeat_structure.repeat_path[0] in unique_flanks or repeat_structure.repeat_path[-1] in unique_flanks:
			print >> tmp_file, 'Alreadu treated this unique flank!'
			print >> tmp_file, repeat_structure
		else:
			repeat_paths[(repeat_structure.repeat_path[0],repeat_structure.repeat_path[-1])] = repeat_structure.repeat_path
			print >> tmp_file, 'Repeat structure taken!'
			print >> tmp_file, repeat_structure
		unique_flanks.add(repeat_structure.repeat_path[0])
		unique_flanks.add(repeat_structure.repeat_path[-1])


		#1 Only take the best repeat path here between unique nodes!!!
		# 2 Also fix so that constraints are not named the same between contig A and 
		# 	contig B since if they are a repeat structure, they can occur multiple times
		# 3 Change the criteria for how we clasify repeats.
		# 4 Lift repeat structure bound threshold??


			# check if repeats can be filled in otherwise, output the repeat structure by itself 
			# as a separate scaffold
			# recalculate positions of contigs
			# make scaffold a la besst1 procedure

	#G = ContigGraph()
	#create_graph(G,contigs)

	print 'LOL NODES left:',len(G.nodes())/2

	besst.get_cut_vertex_cutoff()
	score_file = open(os.path.join(besst.config_params.output_path,'score_info.txt'),'w')
	for lib in besst.libs:
		print >> score_file, lib
	score = 0
	for cut_vertex_size in reversed(metrics.NX_values[-2:]):
		print 'CUTOFF:',cut_vertex_size
		path_factory = paths.PathFactory(besst, G, contigs, cut_vertex_size , 30, 10000)
		#tmp_paths = []
		#print G.edges(data=True)
		for path in path_factory.find_paths():
			print path, path.score, path.good_links,path.bad_links

			for flank1,flank2 in zip(path.path[:-1],path.path[1:]):
				if (flank1,flank2) in repeat_paths:
					#print 'heere__1'
					#print 'repeat path:',repeat_paths[(flank1,flank2)]
					path.add_repeat_region(flank1, flank2, repeat_paths[(flank1,flank2)][1:-1])
					path.find_supporting_links(G_full)
					path.LP_solve_gaps()
					path.score_path(G_full)
					#print 'after:'
					#print path, path.score, path.good_links, path.bad_links, 'lool'

			for flank1,flank2 in zip(reversed(path.path[:-1]),reversed(path.path[1:])):
				if (flank1,flank2) in repeat_paths:
					pass
					#print 'heere__2'
					#print 'repeat path:',repeat_paths[(flank1,flank2)]
					# path.add_repeat_region(flank1, flank2, repeat_paths[(flank1,flank2)][1:-1])
					# path.find_supporting_links(G_full)
					# path.LP_solve_gaps()
					# path.score_path(G_full)
					# print 'after:'
					# print path, path.score, path.good_links, path.bad_links, 'lool'

			# for repeat_ends in repeat_paths:
			# 	if repeat_ends in zip(path.path[:-1],path.path[1:]) or repeat_ends in zip(reversed(path.path[:-1]),reversed(path.path[1:])):
			# 		print 'heere__'
			# 		print 'repeat path:',repeat_paths[repeat_ends]
					# print 'original path :', path
					# path.add_repeat_region(repeat_ends[0],repeat_ends[1],repeat_paths[repeat_ends][1:-1])
					# path.find_supporting_links(G_full)
					# path.LP_solve_gaps()
					# path.score_path(G)
					# print 'after:'
					# print path, path.score, path.good_links, path.bad_links, 'lool'
					# print path.gaps
					# #
			path_factory.add_subpath(path)
			#tmp_paths.append(path)

		scaffolding_score = score_scaffolding(path_factory, G_full)
		print >> score_file, 'CUTOFF:',cut_vertex_size, 'SCORE:', scaffolding_score
		if scaffolding_score > score:
			score = scaffolding_score
			best_scaffolding_path_factory = path_factory

	##
	# Make fasta sequences for all paths
	scaffold_index = 1
	scaffold_file = open(os.path.join(besst.config_params.output_path,'scaffolds.fa'), 'w')
	# for path in tmp_paths:#path_factory.final_paths
	# 	s = sequence.Scaffold(scaffold_index)
	# 	scaffold_index += 1
	# 	scaf = s.make_fasta_string(path,contigs)
	# 	print >> scaffold_file, scaf

	already_printed = set()
	for path in best_scaffolding_path_factory.path_ends.itervalues():#path_factory.final_paths
		end1 = path.path[0][0]
		end2 = path.path[-1][0]
		
		if end1 in already_printed or end2 in already_printed:
			continue
		
		s = sequence.Scaffold('s'+str(scaffold_index))
		scaffold_index += 1
		scaf = s.make_fasta_string(path,contigs,besst.config_params.kmer_overlap)
		print >> scaffold_file, scaf

		already_printed.add(end1)
		already_printed.add(end2)



	for acc, contig_object in contigs.iteritems():
		print acc, contig_object.is_outputted
		if not contig_object.is_outputted: 
			print >> scaffold_file,'>{0}\n{1}'.format(acc,contig_object.sequence)			

	#G.draw(graph_path,repeats)
	#os.subprocess(['dot', '-Tps', '{0}'.format(graph_path),
	#	'{0}'.format(os.path.join(besst.config_params.output_path,'ctg_graph.ps'))])

	## do path searching


	# G = read_input(args)
	# score_list = compute_weighted_interval_solutions(G,args.overlap)
	# make_trusted_paths(G,score_list)
	# # We now expect only linear scaffolds in a perfect world.
	# # We have to remove nodes left with multiple neighbours
	# # TODO: create test instance for this!
	# remove_non_linear_edges(G)
	# scaffold_index = 1
	# scaffold_index = make_scaffolds(G,scaffold_index)

	# #remove all cycles by removing one edge at a time
	# while G.nodes():
	# 	#print 'loop',len(G.nodes())
	# 	remove_cycles(G)
	# 	scaffold_index = make_scaffolds(G,scaffold_index)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="besst-s")
    parser.add_argument('config_file', type=argparse.FileType("r"), help='besst_s config file')
    # parser.add_argument('contigs', type=argparse.FileType("r"), help='contigs fasta file')
    # parser.add_argument('libraries', type=str, nargs='+', help='contigs fasta file')
    # parser.add_argument('output', type=str, help='prefix of output')
    args = parser.parse_args()
    main(args)


