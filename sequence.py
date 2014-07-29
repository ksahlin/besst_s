# class SubSequence(object):
# 	"""Docstring for Subsequence"""
# 	#TODO: Look up what the reverse complement is for all other charachters Y,R,S etc...
# 	reverse_table = {'A':'T','C':'G','G':'C','T':'A', 'N':'N',
# 	'Y':'Y','R':'R','S':'S','M':'M','K':'K','W':'W','n':'n','g':'c','c':'g','t':'a','a':'t'}
# 	def __init__(self, name, contig_object,contig_start,contig_end):
# 		super (SubSequence, self).__init__()
# 		self.name = name
# 		# Contig whereabouts, indexed like lists, i.e.
# 		# 0-indexed and [:7] means positions 0 to 6
# 		self.contig = contig_object
# 		self.contig_start_pos = contig_start
# 		self.contig_end_pos = contig_end

# 		# Scaffold whereabouts, indexed like lists, i.e.
# 		# 0-indexed and [:7] means positions 0 to 6
# 		self.rc = None
# 		self.scaffold = None
# 		self.scaffold_start_pos = None
# 		self.scaffold_end_pos = None
		

# 	def __str__(self):
# 		if self.rc:
# 			return self.reverse_compelment()
# 		else:
# 			return self.contig.sequence[self.contig_start_pos : self.contig_end_pos]

# 	def __repr__(self):
# 		return "From {0}, with coord {1} to {2}".format(self.contig.name,self.contig_start_pos,self.contig_end_pos)

# 	def __len__(self):
# 		return(self.contig_end_pos - self.contig_start_pos)

# 	def __lt__(self,other):
# 		if self.scaffold_start_pos < other.scaffold_start_pos:
# 			return(True)
# 		else:
# 			return(False)

# 	def reverse_complement(self):
# 		string = self.contig.sequence[self.contig_start_pos : self.contig_end_pos]
# 		rc_string = ''.join([SubSequence.reverse_table[nucl] for nucl in reversed(string)])
# 		return(rc_string)


def reverse_complement(string):
 	reverse_table = {'A':'T','C':'G','G':'C','T':'A', 'N':'N',
 	'Y':'Y','R':'R','S':'S','M':'M','K':'K','W':'W','n':'n','g':'c','c':'g','t':'a','a':'t'}
	rc_string = ''.join([reverse_table[nucl] for nucl in reversed(string)])
	return(rc_string)


class Contig(object):
	"""Docstring for Contig"""
	def __init__(self,name, tid=None, sequence = None, length = None):
		super(Contig, self).__init__()
		self.name = name
		self.tid = tid
		self.sequence = sequence
		self.length = len(sequence)
		self.coverage = 0
		self.breakpoints = []
		self.neighbors = {}
		self.subsequences = []
		self.bases_aligned = 0
		self.is_repeat = None 

	def __len__(self):
		return len(self.sequence)

	def get_sequence(self,rc):
		"""
			Returns sequence of contig.
			Reverse complement if rc = 1

			params:
			rc 	bool 
		"""
		if rc == 0:
			return self.sequence
		else:
			return reverse_complement(self.sequence)

	def add_subsequence(self,start,stop,seq):
		self.subsequences.append((start,stop,seq))

	def get_subseq_pos_len(self,contig_pos):
		for start,stop,seq in self.subsequences:
			if contig_pos >= start and contig_pos <= stop:
				return contig_pos - start,seq, stop - start + 1 

		
		return -1,False,-1


class Scaffold(object):
	"""docstring for Scaffold"""
	def __init__(self, name):
		super(Scaffold, self).__init__()
		self.name = name
		#self.subsequences = []

	def make_fasta_string(self,path_object,contigs):
		path = path_object.path
		gap_dict = path_object.gaps
		fasta = '>s{0}\n'.format(self.name)

		# first contig

		fasta += contigs[path[0][0]].get_sequence(1 - path[0][1])

		for ctg1,ctg2 in zip(path[:-1],path[1:]):
			if ctg1[0] == ctg2[0]:
				fasta += contigs[ctg1[0]].get_sequence(ctg1[1])
			else:
				gap = gap_dict[(ctg1,ctg2)]
				print gap
				if gap < 1:
					# do kmer overlap search and merge
					fasta += 'n'
				else:
					fasta += 'N'*int(gap)

		# last contig
		fasta += contigs[ctg2[0]].get_sequence(ctg2[1])

		return fasta

class ScaffoldContainer(object):
	"""docstring for ScaffoldContainer"""
	def __init__(self):
		super(ScaffoldContainer, self).__init__()
		self.scaffolds = {}

	def add_scaffold(self,scaf_object):
		self.scaffolds[scaf_object.name] = scaf_object

	def merge_contig_paths(self,path1,path2):
		"""
			Takes two 
		"""
		
		
