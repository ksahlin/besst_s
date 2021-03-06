import Queue

from pulp import *
from mathstats.normaldist.truncatedskewed import param_est as GC


def other_end(ctg_end):
    if ctg_end[1]:
        return (ctg_end[0],0)
    else:
        return (ctg_end[0],1)

class LinkObservation(object):
    """docstring for LinkObservation"""
    def __init__(self, link_lib, nr_effective_links,mean_obs,from_ctg_index,to_ctg_index, link_type):
        super(LinkObservation, self).__init__()

        self.link_lib           =   link_lib
        self.nr_effective_links =   nr_effective_links
        self.mean_obs           =   mean_obs
        self.from_ctg_index     =   from_ctg_index
        self.to_ctg_index       =   to_ctg_index
        self.link_type          =   link_type
        self.help_variable      =   None

    def gapest_mean_observation(self,contig1,contig2):

        if self.link_type == 'pe':
            self.exp_mean_gapest   =   self.mean_obs + GC.GapEstimator(self.link_lib.mean, self.link_lib.sd, self.link_lib.read_len - self.link_lib.max_softclipped, self.mean_obs, contig1.length, contig2.length)
        if self.link_type == 'pe_am':
            self.exp_mean_gapest   =   self.mean_obs + GC.GapEstimator(self.link_lib.mean_innies, self.link_lib.sd_innies, self.link_lib.read_len - self.link_lib.max_softclipped, self.mean_obs, contig1.length, contig2.length)
        elif self.link_type == 'mp_am' or self.link_type == 'mp':
            self.exp_mean_gapest   =   self.mean_obs + GC.GapEstimator(self.link_lib.mean, self.link_lib.sd, self.link_lib.read_len - self.link_lib.max_softclipped, self.mean_obs, contig1.length, contig2.length)
        
        #print 'mean obs:',self.mean_obs, 'nr_links:',self.nr_effective_links, 'gapest:', self.exp_mean_gapest, 'read_len:',self.link_lib.read_len


class Path(object):
    """docstring for Path"""
    def __init__(self, besst, path, contigs):
        super(Path, self).__init__()
        self.besst = besst
        self.path = path
        self.contigs = contigs
        self.observations = {}
        self.contig_lengths = []

    def __str__(self):
        return str(self.path)

    def __lt__(self,other):
        if self.good_link_ratio < other.good_link_ratio:
            return(True)
        else:
            return(False)

    def score_path(self,graph):
        #self.find_supporting_links(graph)
        #self.LP_solve_gaps()
        try:
            self.good_link_ratio = self.good_links / float(self.bad_links)
        except ZeroDivisionError:
            self.good_link_ratio = self.good_links
        self.score = self.good_link_ratio
        #avg_discrepancy = self.objective / float(self.good_links)
        #self.score = self.good_link_ratio * 1 / (1 + avg_discrepancy)

    def LP_solve_gaps(self):

        for (ctg1,ctg2), link_object in self.observations.iteritems():
            link_object.gapest_mean_observation(self.contigs[ctg1[0]], self.contigs[ctg2[0]])


        #calculate individual exp_mean over each edge given observation with gapest??

        gap_vars= []
        for i in range(len(self.path)/2):
            gap_vars.append( LpVariable(str(i), - self.besst.config_params.kmer_overlap, None, cat='Integer')) #self.mean + 2*self.stddev

        # help variables because objective function is an absolute value

        for (ctg1,ctg2), link_object in self.observations.iteritems():
            link_object.help_variable =  LpVariable("z_"+str(link_object.from_ctg_index)+'_'+str(link_object.to_ctg_index)+'_'+str(link_object.link_type), None, None,cat='Integer')

        problem = LpProblem("PathProblem",LpMinimize)
        problem += lpSum( [ link_object.help_variable*link_object.nr_effective_links for (ctg1,ctg2), link_object in self.observations.iteritems()] ) , "objective"

        # adding constraints induced by the absolute value of objective function
        #print self.observations
        for (ctg1,ctg2), link_object in self.observations.iteritems():
            problem += link_object.exp_mean_gapest - sum(map(lambda x: x, self.contig_lengths[link_object.from_ctg_index+1:link_object.to_ctg_index])) - link_object.mean_obs  - lpSum( gap_vars[ link_object.from_ctg_index : link_object.to_ctg_index] )  <= link_object.help_variable,  "helpcontraint_"+str(link_object.from_ctg_index)+'_'+ str(link_object.to_ctg_index)+'_'+str(link_object.link_type)

        for (ctg1,ctg2), link_object in self.observations.iteritems():
            problem += - link_object.exp_mean_gapest + sum(map(lambda x: x, self.contig_lengths[link_object.from_ctg_index+1:link_object.to_ctg_index])) + link_object.mean_obs + lpSum( gap_vars[ link_object.from_ctg_index : link_object.to_ctg_index] )  <= link_object.help_variable,  "helpcontraint_negative_"+str(link_object.from_ctg_index)+'_'+ str(link_object.to_ctg_index)+'_'+str(link_object.link_type)

        try:
            problem.solve()
        except : #PulpSolverError:
            problem.solve()
            print 'Could not solve LP, printing instance:'
            print 'Objective:'
            print problem.objective
            print 'Constraints:'
            print problem.constraints


        optimal_gap_solution = [0]*( len(self.path)/2)
        self.gaps = {}
        #print self.path
        for v in problem.variables():
            try:
                self.gaps[(self.path[int(v.name)*2], self.path[int(v.name)*2+1]) ] = v.varValue #self.path.index(ctg)/2 , self.path.index(nbr)/2 + 1
                self.gaps[(self.path[int(v.name)*2+1], self.path[int(v.name)*2]) ] = v.varValue
                optimal_gap_solution[int( v.name)] = v.varValue
                #print v.name, "=", v.varValue
            except ValueError:
                #print v.name, "=", v.varValue,
                pass

        self.objective = value(problem.objective)
        #print 'hej:',self.gaps #self.objective,optimal_gap_solution
        return optimal_gap_solution       

    def find_supporting_links(self,graph):
        self.good_links = 0 
        self.bad_links = 0
        path_forward, path_backwards = [], []
        for i, ctg in enumerate(self.path):
            path_forward.append(ctg) if i %2 ==0 else path_backwards.append(ctg)

        for ctg in path_forward:
            #print 'ctg:',ctg
            for nbr in graph.neighbors(ctg):
                #print 'nbr:',nbr
                if nbr not in path_backwards or self.path.index(nbr) < self.path.index(ctg):
                    for edge in graph.all_edges_between_two_nodes(ctg,nbr):
                        self.bad_links += edge['n']
                        
                else:
                    for edge in graph.all_edges_between_two_nodes(ctg,nbr):
                        self.good_links += edge['n']
                        self.observations[(ctg,nbr)] = LinkObservation(edge['l'], graph.effective_nr_of_links(edge),edge['o_bar'],self.path.index(ctg)/2 , self.path.index(nbr)/2 + 1, edge['t'])
                        self.contig_lengths.append(self.contigs[ctg[0]].length)

        #adding last contig length
        self.contig_lengths.append(self.contigs[self.path[-1][0]].length)

        for ctg in path_backwards:
            for nbr in graph.neighbors(ctg):
                if nbr not in path_forward:
                    for edge in graph.all_edges_between_two_nodes(ctg,nbr):
                        self.bad_links += edge['n']
                else:
                    pass


    def add_repeat_region(self, flank1,flank2,repeats):
        pos1 = self.path.index(flank1)
        pos2 = self.path.index(flank2)
        assert abs(pos1-pos2) == 1
        if pos1 < pos2:
            self.path[pos1+1:pos1+1] = repeats
        else:
            self.path[pos2+1:pos2+1] = repeats[::-1]


class PathFactory(object):
    """docstring for PathFactory"""
    def __init__(self, besst, graph, contigs, cut_length, max_depth, iter_tresh):
        super(PathFactory, self).__init__()
        self.besst = besst
        self.graph = graph
        self.contigs = contigs
        self.max_depth = max_depth
        self.iter_tresh = iter_tresh
        self.cut_vertices = set()

        self.path_ends = {}
        self.already_merged = set()

        for ctg in filter(lambda ctg: contigs[ctg].length > cut_length and (contigs[ctg].name,0) in self.graph, contigs):
            self.cut_vertices.add((contigs[ctg].name, 0))
            self.cut_vertices.add((contigs[ctg].name,1))
        self.paths = []


    def find_paths(self):
        self.already_visited = set()
        for start in self.cut_vertices:
            high_score_paths = []
            self.forbidden = set([self.graph.other_end(start)])
            print 'treating contig:', start
            contig_paths = []
            for path in self.BFS(start):
                p = Path(self.besst, path, self.contigs)
                #p.score_path(self.graph)
                #self.paths.append(p)
                p.find_supporting_links(self.graph)
                p.score_path(self.graph)
                contig_paths.append(p)
            sorted_paths = sorted(contig_paths, reverse=True)
            try:
                high_score_paths = sorted_paths[0:5]
            except:
                continue

            #print high_score_path.good_links,high_score_path.bad_links
            for path in high_score_paths:
                path.LP_solve_gaps()
                self.paths.append(path) 
            self.already_visited.add(start)
        return sorted(self.paths, reverse=True)

    def BFS(self,start):
        
        temp_path = [start]
        self.queue = Queue.Queue() #MyQUEUE()
        self.queue.put(temp_path)

        i = 0
        while not self.queue.empty() and i <= self.iter_tresh:
            #i+=1
            #if i % 5000 == 0:
            #    print 'BFS:',i
            tmp_path = self.queue.get() 
            last_node = tmp_path[-1]
            try:
                second_last_node = tmp_path[-2]
            except IndexError:
                second_last_node = tmp_path[-1]

            if len(tmp_path) > self.max_depth or last_node in self.forbidden or second_last_node in self.already_visited:
                continue

            if last_node != start and last_node in self.cut_vertices:
                yield tmp_path[:-1]
                continue


            for link_node in self.graph.neighbors(last_node):
                if link_node not in tmp_path and self.graph.other_end(link_node) not in tmp_path:
                    new_path = []
                    new_path = tmp_path + [link_node, self.graph.other_end(link_node)]
                    self.queue.put(new_path)
   

    def merge_contig_paths(self,path1,path2):
        """
            Takes two path path1 and path2 and merges them
        """

        ##
        # dectect in what ends they should be merged

        ## rare circular case
        if path1.path[0] == other_end(path2.path[0]) and  path1.path[-1] == other_end(path2.path[-1]) \
        or path1.path[-1] == other_end(path2.path[0]) and  path1.path[0] == other_end(path2.path[-1]):
            print ' circular paths! Skipping merging'
            return path1

        elif path1.path[-1] == other_end(path2.path[0]) :
            p = Path(self.besst,path1.path + path2.path, self.contigs )
            p.gaps = dict(path1.gaps.items() + path2.gaps.items())
            return p
               
        elif path1.path[-1] == other_end(path2.path[-1]):
            p = Path(self.besst,path1.path + [i for i in reversed(path2.path)], self.contigs )
            p.gaps = dict(path1.gaps.items() + path2.gaps.items())
            return p       


        elif path1.path[0] == other_end(path2.path[0]):
            p = Path(self.besst, [i for i in reversed(path1.path)] + path2.path, self.contigs )
            p.gaps = dict(path1.gaps.items() + path2.gaps.items())
            return p            

        elif path1.path[0] == other_end(path2.path[-1]):
            p = Path(self.besst, path2.path + path1.path, self.contigs )
            p.gaps = dict(path1.gaps.items() + path2.gaps.items()) 
            return p           

        else:

            print 'NOOOO'
            return 0


    def add_subpath(self,path):

        ##
        # There exists a unique (internal) contig that is included in another path already

        for node in path.path:
            if not self.contigs[node[0]].is_repeat and node in self.already_merged:
                #print 'Unique ctg already treated'
                return



        ##
        # look if endpoint in current path is the other flank of a contig that is 
        # a flank of another path that has been treated previously. If this is the case,
        # merge.

        if other_end(path.path[0]) in self.path_ends and other_end(path.path[-1]) in self.path_ends:
             new_path_temp = self.merge_contig_paths(path, self.path_ends[other_end(path.path[0])])
             new_path = self.merge_contig_paths(new_path_temp, self.path_ends[other_end(path.path[-1])])
             del self.path_ends[other_end(path.path[0])]
             del self.path_ends[other_end(path.path[-1])]

        elif other_end(path.path[0]) in self.path_ends:
             new_path = self.merge_contig_paths(path, self.path_ends[other_end(path.path[0])])
             del self.path_ends[other_end(path.path[0])]

        elif other_end(path.path[-1]) in self.path_ends:
            new_path = self.merge_contig_paths(path, self.path_ends[other_end(path.path[-1])])
            del self.path_ends[other_end(path.path[-1])]

        else:
            new_path = path

        ##
        # Add all treated unique internal contigs to already visited
        for node in new_path.path[1:-1]:
            if not self.contigs[node[0]].is_repeat:
                self.already_merged.add(node)     

        ##
        # Add the endpoints of new_path path to dictionary of path ends.
        self.path_ends[ new_path.path[0] ] = new_path
        self.path_ends[ new_path.path[-1] ] = new_path


# class MyQUEUE(object): # just an implementation of a queue
    
#     def __init__(self):
#         self.holder = []
        
#     def enqueue(self,val):
#         self.holder.append(val)
        
#     def dequeue(self):
#         val = None
#         try:
#             val = self.holder[0]
#             if len(self.holder) == 1:
#                 self.holder = []
#             else:
#                 self.holder = self.holder[1:]   
#         except:
#             pass
            
#         return val  
        
#     def IsEmpty(self):
#         result = False
#         if len(self.holder) == 0:
#             result = True
#         return result




