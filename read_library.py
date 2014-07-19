import sys

from mathstats.normaldist.normal import MaxObsDistr

import bam_parser


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
	def __init__(self, name, bam_path, lib_type):
		super(Library, self).__init__()
		self.name = name
		self.bam_path = bam_path
		self.lib_type = lib_type


	def get_libary_metrics(self):
		ins_size_reads_innies = []
		ins_size_reads_outies = []
		for read_type,read in bam_parser.BamParser(self.bam_path).proper_aligned_unique_pairs('bwa_mem',samples=1000000):
			if self.lib_type == 'MP': 
				if read_type == 'innie':
					ins_size_reads_innies.append(abs(read.tlen))
				elif read_type == 'outie':
					ins_size_reads_outies.append(abs(read.tlen))
			#print read.tlen #, read
		#self.mean =
		#self.sd =
		self.mean_innies, self.sd_innies = remove_outliers(ins_size_reads_innies)
		self.mean_outies, self.sd_outies = remove_outliers(ins_size_reads_outies)

		count_contamine, count_total = 0,0

		if self.lib_type == 'MP':
			for read in bam_parser.BamParser(self.bam_path).aligned_reads():
				if bam_parser.is_proper_aligned_unique_innie(read):
					count_contamine += 1
				count_total += 1
			self.contamine_rate = count_contamine/float(count_total)
			print self.contamine_rate
		else:
			self.contamine_rate = False





def main(bam_path):
	lib = Library('test',bam_path,'MP')
	lib.get_libary_metrics()

if __name__ == '__main__':
	main(sys.argv[1])