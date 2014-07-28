def get_contigs(contig_file):
    k = 0
    sequence = ''
    accession = ''
    for line in contig_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            k += 1
        elif line[0] == '>':
            if len(sequence) > 0:
        	   yield accession, sequence
            sequence = ''
            accession = line[1:].strip().split()[0]
        else:
            sequence += line.strip()
    
    if len(sequence) > 0:
        yield accession, sequence