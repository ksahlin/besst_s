import os, sys
from operator import itemgetter
from collections import defaultdict

from mathstats.normaldist.normal import MaxObsDistr

from parser import bam_parser


def AdjustInsertsizeDist(mean_insert, std_dev_insert, insert_list):
    k = MaxObsDistr(len(insert_list), 0.95)
    filtered_list = list(filter((lambda x : (x < mean_insert + k * std_dev_insert and x > mean_insert - k * std_dev_insert)), insert_list))
    if len(insert_list) > len(filtered_list):
        return(True, filtered_list)
    else:
        return(False, filtered_list)

def remove_outliers(ins_size_reads):
	## SMOOTH OUT THE MEAN HERE by removing extreme observations
	if not ins_size_reads:
		return 0,0
	n = float(len(ins_size_reads))
	mean_isize = sum(ins_size_reads) / n
	std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5
	extreme_obs_occur = True
	while extreme_obs_occur:
	    extreme_obs_occur, filtered_list = AdjustInsertsizeDist(mean_isize, std_dev_isize, ins_size_reads)
	    n = float(len(filtered_list))
	    mean_isize = sum(filtered_list) / n
	    std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), filtered_list))) / (n - 1)) ** 0.5
	    ins_size_reads = filtered_list
	n = float(len(ins_size_reads))
	mean_isize = sum(ins_size_reads) / n
	std_dev_isize = (sum(list(map((lambda x: x ** 2 - 2 * x * mean_isize + mean_isize ** 2), ins_size_reads))) / (n - 1)) ** 0.5
	return mean_isize, std_dev_isize


class Library(object):
	"""docstring for Library"""
	def __init__(self, lib_name, lib_type, aligner, bam_path, link_file_path, bam_path2=None, mean=None,sd=None):
		super(Library, self).__init__()
		self.lib_name = lib_name
		self.bam_path = bam_path
		self.bam_path2 = bam_path2
		self.lib_type = lib_type
		self.aligner = aligner
		#open(link_file_path,'w').close()
		#os.utime(link_file_path, None)
		self.link_file_path = link_file_path
		self.mean=mean
		self.sd=sd


	def calculate_coverage(self, contigs):
		bam_iter = bam_parser.BamParser(self.bam_path,bam_path2=self.bam_path2)

		self.max_softclipped = 0

		if self.aligner == 'bowtie':
			for read1,read2 in bam_iter.aligned_reads(self.aligner):
				if read1:
					contigs[bam_iter.bam_file.getrname(read1.tid)].bases_aligned += read1.alen
					if read1.rlen - read1.qlen > self.max_softclipped:
						self.max_softclipped = read1.rlen - read1.qlen

				if read2:
					contigs[bam_iter.bam_file.getrname(read2.tid)].bases_aligned += read2.alen
					if read2.rlen - read2.qlen > self.max_softclipped:
						self.max_softclipped = read2.rlen - read2.qlen

		elif self.aligner == 'bwamem' or self.aligner == 'bwa':
			for read in bam_iter.aligned_reads(self.aligner):
				contigs[bam_iter.bam_file.getrname(read.tid)].bases_aligned += read.alen
				if read.rlen - read.qlen > self.max_softclipped:
					self.max_softclipped = read.rlen - read.qlen

		lib_cov = []
		for contig in contigs.itervalues():
			contig.coverage += contig.bases_aligned / float(contig.length)
			lib_cov.append(contig.bases_aligned / float(contig.length))
			contig.bases_aligned = 0
		
		self.mean_coverage = sum([cov for cov in lib_cov]) / len(contigs)
		self.sd_coverage = ( sum([ (cov - self.mean_coverage)**2 for cov in lib_cov]) / len(contigs) )**0.5

	def get_libary_metrics(self):
		ins_size_reads_innies = []
		ins_size_reads_outies = []
		read_length = []
		for read_type,read in bam_parser.BamParser(self.bam_path,bam_path2=self.bam_path2).proper_aligned_unique_pairs(self.aligner,samples=1000000):
			if self.lib_type == 'mp': 
				if read_type == 'innie':
					ins_size_reads_innies.append(abs(read.tlen))
				elif read_type == 'outie':
					ins_size_reads_outies.append(abs(read.tlen))
			elif self.lib_type =='pe':
				ins_size_reads_innies.append(abs(read.tlen))

			read_length.append(read.qlen)

		self.mean_innies, self.sd_innies = remove_outliers(ins_size_reads_innies)
		self.mean_outies, self.sd_outies = remove_outliers(ins_size_reads_outies)

		count_contamine, count_total = 0,0

		if self.lib_type == 'mp' and self.aligner == 'bowtie':
			for read1,read2 in bam_parser.BamParser(self.bam_path,bam_path2=self.bam_path2).aligned_reads(self.aligner):
				if read1 and read2 and bam_parser.proper_unique_alignment_innie_bowtie(read1,read2):
					count_contamine += 1
				elif read1 and read2:
					count_total += 2
				elif read1:
					count_total += 1
				elif read2:
					count_total += 1
			self.contamine_rate = count_contamine/float(count_total)
		

		elif self.lib_type == 'mp' and (self.aligner == 'bwa' or self.aligner == 'bwamem'):
			for read in bam_parser.BamParser(self.bam_path,bam_path2=self.bam_path2).aligned_reads(self.aligner):
				if bam_parser.is_proper_aligned_unique_innie(read):
					count_contamine += 1
				count_total += 1
			self.contamine_rate = count_contamine/float(count_total)
		else:
			self.contamine_rate = False

		self.read_len = sum(read_length)/len(read_length)

		if not self.mean and not self.sd:
			self.mean = self.mean_innies if self.lib_type =='pe' else self.mean_outies
			self.sd = self.sd_innies if self.lib_type =='pe' else self.sd_outies

	def sorted_bam_to_link_file(self,contigs, kmer_overlap):
		bam_iter = bam_parser.BamParser(self.bam_path,bam_path2=self.bam_path2)
		edges = {}
		for contig in contigs.keys():
			edges[contig] = ContigConnections(contig)

		for read1,read2 in bam_iter.unique_reads_on_different_references(self.aligner):
			if (self.aligner == 'bwamem' or self.aligner == 'bwa') and (read1.rnext != read2.tid or read2.rnext != read1.tid):
				print read1,read2
			#print read1,read2
			index, ctg = min(enumerate([read1.tid,read2.tid]), key=itemgetter(1))
			#print bam_iter.bam_file.getrname(read1.tid),bam_iter.bam_file.getrname(read2.tid)
			if self.lib_type == 'mp':
				mp_obs1, mp_obs2 = bam_parser.get_mp_observation(read1, read2, bam_iter.contig_lengths[ bam_iter.bam_file.getrname(read1.tid)],bam_iter.contig_lengths[ bam_iter.bam_file.getrname(read2.tid)] )
				pe_obs1, pe_obs2 = bam_parser.get_pe_observation(read1, read2, bam_iter.contig_lengths[ bam_iter.bam_file.getrname(read1.tid)],bam_iter.contig_lengths[ bam_iter.bam_file.getrname(read2.tid)] )
				
				##
				# reads could be both RF or FR, we can't deduce because of the size of the contigs  
				if mp_obs1 + mp_obs2 < self.mean + 4*self.sd + kmer_overlap and pe_obs1 + pe_obs2 < self.mean_innies + 4*self.sd_innies + kmer_overlap:
					if index == 0:
						edges[bam_iter.bam_file.getrname(ctg)].add_pe_link(pe_obs1,pe_obs2, read1.is_reverse, read2.is_reverse, bam_iter.bam_file.getrname(read2.tid),'pe_am') 
						edges[bam_iter.bam_file.getrname(ctg)].add_mp_link(mp_obs1,mp_obs2, read1.is_reverse, read2.is_reverse, bam_iter.bam_file.getrname(read2.tid),'mp_am')
					else:
						edges[bam_iter.bam_file.getrname(ctg)].add_pe_link(pe_obs2,pe_obs1, read2.is_reverse, read1.is_reverse, bam_iter.bam_file.getrname(read1.tid),'pe_am') 
						edges[bam_iter.bam_file.getrname(ctg)].add_mp_link(mp_obs2,mp_obs1, read2.is_reverse, read1.is_reverse, bam_iter.bam_file.getrname(read1.tid),'mp_am')
			
				##
				# Reads can only be RF
				elif mp_obs1 + mp_obs2 < self.mean + 4*self.sd + kmer_overlap :	
					if index == 0:
						edges[bam_iter.bam_file.getrname(ctg)].add_mp_link(mp_obs1,mp_obs2, read1.is_reverse, read2.is_reverse, bam_iter.bam_file.getrname(read2.tid), 'mp') 
					else:
						edges[bam_iter.bam_file.getrname(ctg)].add_mp_link(mp_obs2,mp_obs1, read2.is_reverse, read1.is_reverse, bam_iter.bam_file.getrname(read1.tid),'mp') 
					
				##
				# Reds must be FR (contamine reads)

				elif pe_obs1 + pe_obs2 < self.mean_innies + 4*self.sd_innies + kmer_overlap:
					if index == 0:
						edges[bam_iter.bam_file.getrname(ctg)].add_pe_link(pe_obs1,pe_obs2, read1.is_reverse, read2.is_reverse, bam_iter.bam_file.getrname(read2.tid),'pe') 
					else:
						edges[bam_iter.bam_file.getrname(ctg)].add_pe_link(pe_obs2,pe_obs1, read2.is_reverse, read1.is_reverse, bam_iter.bam_file.getrname(read1.tid), 'pe') 
				else:
					print 'spurious link'

			else:
				obs1, obs2 = bam_parser.get_pe_observation(read1, read2, bam_iter.contig_lengths[ bam_iter.bam_file.getrname(read1.tid)], bam_iter.contig_lengths[ bam_iter.bam_file.getrname(read2.tid)])
				if obs1 + obs2 < self.mean + 4*self.sd + kmer_overlap:
					if index == 0:
						edges[bam_iter.bam_file.getrname(ctg)].add_pe_link(obs1,obs2, read1.is_reverse, read2.is_reverse, bam_iter.bam_file.getrname(read2.tid),'pe') 
					else:
						edges[bam_iter.bam_file.getrname(ctg)].add_pe_link(obs2, obs1, read2.is_reverse, read1.is_reverse, bam_iter.bam_file.getrname(read1.tid),'pe') 

			 #,  bam_iter.contig_lengths[ bam_iter.bam_file.getrname(read1.tid)],bam_iter.contig_lengths[ bam_iter.bam_file.getrname(read2.tid)]
			#print read1,read2 #, read.tlen
			#link_file_path.write(str(seq_id)+'\t'+contig.name+'\t'+ str(start)+'\t'+str(stop)+'\n')

		link_file = open(self.link_file_path,'w')
		for connection in edges:
			if str(edges[connection]):
				link_file.write(str(edges[connection]))

		link_file.close()
		#self.link_file.seek(0)
		#self.link_file.read()
				#print >> self.link_file, edges[connection]
	

	def __str__(self):
		if self.lib_type == 'mp':
			return 'Library type:{0}\nLibrary mean:{1} \
			\nLibrary sd:{2}\nAligned with:{3}\nContamination rate:{4}\nContamination mean:{5}\n Contamination sd:{6}\nLibrary mean cov:{7}\nLibrary sd cov:{8}'.format(self.lib_type, 
				self.mean,self.sd,self.aligner,
				self.contamine_rate, self.mean_innies,self.sd_innies,
				self.mean_coverage,self.sd_coverage)
		else:
			return 'Library type:{0}\nLibrary mean:{1}\nLibrary sd:{2}\nAligned with:{3}\nLibrary mean cov:{4}\nLibrary sd cov:{5}'.format(self.lib_type, self.mean,self.sd,self.aligner,self.mean_coverage,self.sd_coverage)



class ContigConnections(object):
    """docstring for ContigConnections"""
    def __init__(self, name):
        super(ContigConnections, self).__init__()
        self.name = name
        self.connections_forward = defaultdict(list)
        self.connections_reverse = defaultdict(list)     

    def add_mp_link(self, obs1, obs2, is_rc1, is_rc2, other,link_type): 
        if is_rc1:
            self.connections_forward[(other,int(is_rc2),link_type)].append((obs1,obs2))
        else:
            self.connections_reverse[(other,int(is_rc2), link_type)].append((obs1,obs2))
    def add_pe_link(self, obs1, obs2, is_rc1, is_rc2,other,link_type ):
        if not is_rc1:
            self.connections_forward[(other,1-is_rc2,link_type)].append((obs1,obs2))
        else:
            self.connections_reverse[(other,1-is_rc2, link_type)].append((obs1,obs2))


    def make_connections_to_string(self,dict_,not_rc):
        i_am = ''
        for edge in dict_:
			#print self.name, not_rc, edge[0], str(edge[1]), len(dict_[edge]),str(edge[2])
			i_am += '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(self.name, not_rc, edge[0], str(edge[1]), len(dict_[edge]),str(edge[2])) #replace 0 with gap
			obs_first = '# '
			obs_second = '# '
			for link in dict_[edge]:
			    obs_first += str(link[0]) + ' '
			    obs_second += str(link[1]) + ' '
			    
			i_am += obs_first+'\n' + obs_second +'\n'
        return(i_am)

    def __str__(self):
    	#print 'final tables: ', self.connections_forward, self.connections_reverse
    	#print 'forward:', self.connections_forward
    	#print 'reverse:', self.connections_reverse
        return self.make_connections_to_string(self.connections_forward,1) + self.make_connections_to_string(self.connections_reverse,0)



def main(bam_path):
	lib = Library('test',bam_path,'MP')
	lib.get_libary_metrics()

if __name__ == '__main__':
	main(sys.argv[1])