import itertools

def get_orientation(o,s1,s2):
    if int(o) == 1:
        return True 
    elif int(o) == 0:
        return False
    else:
        sys.stderr.write('Orientation incorrectly specified for link: {0}, {1}.'.format(s1,s2))
        raise IOError

def get_links(link_file_path):
    link_file = open(link_file_path,'r')
    while True:
        next_3_lines = list(itertools.islice(link_file, 3))
        if not next_3_lines:
            break
        else:
            ctg1, o1, ctg2, o2, nr_links, link_type = next_3_lines[0].split()
            obs_list1 = map(lambda x: int(x), next_3_lines[1][1:].split())
            obs_list2 = map(lambda x: int(x), next_3_lines[2][1:].split())  
            mean_obs = sum(map(lambda x: x[0]+x[1], zip(obs_list1,obs_list2)) )/float(nr_links)  
            yield (ctg1, int(o1), ctg2, int(o2), int(nr_links), link_type, mean_obs)

    #link_file.seek(0)
    link_file.close()