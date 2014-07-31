import os
import sys

class ConfigParams(object):
    """docstring for ConfigParams"""
    def __init__(self):
        super(ConfigParams, self).__init__()
        self.libs = []

    def read_cfg(self,config_file):
        for line in config_file:
            if line[:1] == "1":
                self.output_path = line.split()[1] 
                if not os.path.exists(self.output_path):
   	 				os.makedirs(self.output_path)
            
            elif line[:1] == "2":
            	if not os.path.exists(line.split()[1]):
            		sys.stderr.write('Contig path given in config file does not exist. You gave:\n {0}\n'.format(line.split()[1]))
            		sys.exit()
                self.contig_file ,self.kmer_overlap = open(line.split()[1],'r'), int(line.split()[2])

            
            elif line[:1] == "3":
                input_lib = line.split()
                lib_type, aligner, lib_location  =  input_lib[1], input_lib[2], input_lib[3]
            	if not os.path.exists(lib_location):
            		sys.stderr.write('Library path given in config file does not exist. You gave:\n {0}\n'.format(lib_location))
            		sys.exit()
                if len(input_lib) == 6:
                    mean,sd = int(input_lib[-2]), int(input_lib[-1])
                    lib = (lib_type, aligner, lib_location, mean, sd)
                elif len(input_lib) == 4:
                    lib = (lib_type, aligner, lib_location)
                self.libs.append(lib)

            elif line[:1] == "4":
            	self.min_links = int(line.split()[1])
