
class AssemblyMetrics(object):
	"""docstring for AssemblyMetrics"""
	def __init__(self, sequences):
		super(AssemblyMetrics, self).__init__()
		self.sequences = sequences
	
	def __str__(self):
		string = ''
		for x in range(10,101,10):
			string += 'N{0} = {1}\n'.format(x,self.NX(x))
		for x in range(10,101,10):
			string += 'L{0} = {1}\n'.format(x,self.LX(x))

		string += 'E-size: {0}\n'.format(self.E_size())
		string += 'Longest contig: {0}\n'.format(self.longest_contig())
		string += 'Ratio of N\'s in assembly: {0}\n'.format(self.N_ratio())
		string += 'Total assembly length: {0}\n'.format(self.total_length())
		return string

	def LX(self,x):
		lengths = map(lambda x: len(self.sequences[x]), self.sequences )
		lengths.sort(reverse=True)
		L100 = sum(lengths)
		stop = L100/ (100/float(x))
		curr_sum = 0
		nr_ctgs = 0
		for length in lengths:
			curr_sum += length
			nr_ctgs += 1
			if curr_sum >= stop:
				return nr_ctgs
	def NX(self,x):
		lengths = map(lambda x: len(self.sequences[x]), self.sequences )
		lengths.sort(reverse=True)
		N100 = sum(lengths)
		stop = N100/ (100/float(x))
		curr_sum = 0
		for length in lengths:
			curr_sum += length
			if curr_sum >= stop:
				return length

	def E_size(self):
		lengths = map(lambda x: len(self.sequences[x]), self.sequences )
		G = sum(lengths) 
		contig_square_sum = sum(map(lambda x: x**2, lengths))
		return contig_square_sum /float(G)

	def longest_contig(self):
		lengths = map(lambda x: len(self.sequences[x]), self.sequences )
		lengths.sort(reverse=True)
		return lengths[0]

	def N_ratio(self):
		return sum( map(lambda x: self.sequences[x].sequence.count('N'), self.sequences ))

	def total_length(self):
		return sum(map(lambda x: len(self.sequences[x]), self.sequences ))


	def coverage(self):
		return