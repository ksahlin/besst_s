import argparse
import sys
import os

from parser import fasta_parser,config_parser,link_parser
import sequence
import read_library
from contig_graph import ContigGraph 

class GlobaInfo(object):
	"""docstring for GlobaInfo"""
	def __init__(self):
		super(GlobaInfo, self).__init__()
		self.config_params = None
		self.libs = []

	def __str__(self):
		return

besst = GlobaInfo()

def read_in_contigs(config_params,contig_dict):
	for accession,seq in fasta_parser.get_contigs(config_params.contig_file):
		c = sequence.Contig(accession, sequence = seq)
		contig_dict[c.name] = c

def read_library_metrics(config_params, contigs):
	lib_name =1
	for lib_type, aligner, lib_loc in config_params.libs:
		lib = read_library.Library(lib_type, aligner, lib_loc, os.path.join(config_params.output_path, 'lib_{0}_links.txt'.format(lib_name)))
		lib.get_libary_metrics()
		print contigs
		lib.calculate_coverage(contigs)
		lib_name += 1
		besst.libs.append(lib)
		print lib

		

def create_link_file(contigs):
	for lib in besst.libs:
		lib.sorted_bam_to_link_file(contigs,besst.config_params.kmer_overlap)


def create_graph():
	G = ContigGraph()


	##
	# Create edges
	for lib in besst.libs:
		for (ctg1, orientation1, ctg2, orientation2, link_count, link_type, mean_obs) in  link_parser.get_links(lib.link_file):
			G.add_edge(ctg1, orientation1, ctg2, orientation2, link_count, link_type, mean_obs, lib)
	print G
				#G.add_edge((subsequences[ seq1 ] ,orientation1),(subsequences[ seq2 ],orientation2),d=naive_gap,s=score.nr_links(link_count))

	G.remove_self_links()

# 	return(G)


def main(args):
	params = config_parser.ConfigParams()
	params.read_cfg(args.config_file)
	besst.config_params = params
	contigs = {}
	read_in_contigs(params,contigs)
	read_library_metrics(params,contigs)
	create_link_file(contigs)
	create_graph()
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


