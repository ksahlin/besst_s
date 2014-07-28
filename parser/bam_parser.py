import os
import sys
import pysam

##
# Opens a .bam or .sam file and returns the file
# object.
#
# @param bam_file_path Path to the .bam or .sam.
# 
# @return File object for .bam or .sam.
#
def open_bam_file(bam_file_path):
    bam_file_name, bam_file_ext = os.path.splitext(bam_file_path)
    if bam_file_ext == ".bam":
        return pysam.Samfile(bam_file_path, 'rb')
    elif bam_file_ext == ".sam":
        return pysam.Samfile(bam_file_path, 'r')
    else:
        return IOError("open_bam_file: File must be either .bam or .sam.")

def is_proper_aligned_unique_innie(read):
    return (read.is_reverse and not read.mate_is_reverse and  read.tlen < 0 and read.rname == read.mrnm) or \
                (not read.is_reverse and read.mate_is_reverse and read.is_read2 and read.tlen > 0 and read.rname == read.mrnm ) \
                and not read.mate_is_unmapped and read.opt('XT')=='U' and not read.is_secondary
def is_proper_aligned_unique_outie(read):
    return (read.is_reverse and not read.mate_is_reverse and  read.tlen > 0 and read.rname == read.mrnm) or \
                (not read.is_reverse and read.mate_is_reverse and read.is_read2 and read.tlen < 0 and read.rname == read.mrnm ) \
                and not read.mate_is_unmapped and read.opt('XT')=='U' and not read.is_secondary

def is_unique_read_link(read):
    return not read.is_unmapped and not read.mate_is_unmapped and read.rname != read.mrnm \
    and read.opt('XT')=='U' and not read.is_secondary


def get_mp_observation(read1, read2, ctg_len1, ctg_len2):
    if not read1.is_reverse:
        obs1 = read1.pos + read1.qlen
    else:
        obs1 = ctg_len1 - read1.pos

    if not read2.is_reverse:
        obs2 = read2.pos + read2.qlen
    else:
        obs2 = ctg_len2 - read2.pos
    return obs1, obs2

def get_pe_observation(read1, read2, ctg_len1, ctg_len2):
    if read1.is_reverse:
        obs1 = read1.pos + read1.qlen
    else:
        obs1 = ctg_len1 - read1.pos

    if read2.is_reverse:
        obs2 = read2.pos + read2.qlen
    else:
        obs2 = ctg_len2 - read2.pos
    return obs1, obs2

class BamParser(object):
    """docstring for BamParser"""
    def __init__(self, bam_file):
        super(BamParser, self).__init__()
        self.bam_file = open_bam_file(bam_file)
        self.contig_lengths = dict(zip(self.bam_file.references,self.bam_file.lengths))

    def proper_aligned_unique_pairs(self,aligner, samples=2**32):
        nr_samples = 0
        if aligner == 'bwa' or aligner == 'bwa_mem':
            for ref, length in sorted(zip(self.bam_file.references,self.bam_file.lengths),key = lambda x: x[1], reverse = True):
                try:
                    iter_ = self.bam_file.fetch(ref)
                except ValueError:
                    sys.stderr.write('Need indexed bamfiles, index file should be located in the same directory as the BAM file\nterminating..\n')
                    sys.exit(0)
                for read in iter_:
                    if read.is_read2 and is_proper_aligned_unique_innie(read):
                        nr_samples += 1
                        yield 'innie',read
                    elif read.is_read2 and is_proper_aligned_unique_outie(read): 
                        nr_samples += 1
                        yield 'outie',read

                    if nr_samples >= samples:
                        break
        self.bam_file.seek(0)

    def aligned_reads(self):
        for read in self.bam_file:
            if not read.is_unmapped:
                yield read 
        self.bam_file.seek(0)

    def unique_reads_on_different_references(self, aligner):
        read_pairs = {}
        if aligner == 'bwa_mem':
            for read in self.bam_file:
                if is_unique_read_link(read):
                    #tmp for tests:
                    #print read.qname[:-1]
                    #print read_pairs
                    if read.qname[:-1] in read_pairs:
                        #print 'lol'
                        #read2 = read_pairs[read.qname]
                        #if read.tid == read2.tid:    
                        yield read, read_pairs[read.qname[:-1]]
                        #else: 
                        #    pass
                        del read_pairs[read.qname[:-1]]
                    else:
                        read_pairs[read.qname[:-1]] = read
                    


        elif aligner == 'bwa':
            pass
        elif aligner == 'bowtie':
            pass
        self.bam_file.seek(0)

    def long_reads_for_coverage(self):
        pass
    def long_reads_scaffold_links(self):
        pass

def main(path):
    for read in BamParser(path).reads_for_coverage():
        if read.qlen != 100:
        #print read.tags
            print read.qlen
        #print read.opt('XT')


if __name__ == '__main__':
    main(sys.argv[1])
