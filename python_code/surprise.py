#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
computeSurprise calculates the Surprise of the partition of a network
Surprise parameters are F, M, n and p

Copyright (C) 2012 Rodrigo Aldecoa and Ignacio Marín                   
                                                                        
 This program is free software: you can redistribute it and/or modify   
 it under the terms of the GNU General Public License as published by   
 the Free Software Foundation, either version 3 of the License, or      
 (at your option) any later version.                                    
                                                                        
 This program is distributed in the hope that it will be useful,        
 but WITHOUT ANY WARRANTY; without even the implied warranty of         
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 GNU General Public License for more details.                           
                                                                        
 You should have received a copy of the GNU General Public License      
 along with this program.  If not, see <http://www.gnu.org/licenses/>.  
                                                                        
Contact info: Rodrigo Aldecoa <raldecoa@ibv.csic.es>                    
Converted to python-igraph by @datapornstar <datapornstar@riseup.net> 

 If you use this program, please cite:
       Aldecoa R, Marín I(2011)
       Deciphering network community structure by Surprise
       PLoS ONE 6(9): e24195

"""

from sys import argv 
from math import log, pow as mpow

def compute_surprise(F, M, n, p):
    j = p
    logP = log_hyper_probability(F, M, n, j)

    minimum = M
    if n < M: minimum = n

    is_enough = False
    while not is_enough and (j < minimum):
        j = j + 1
        next_logP = log_hyper_probability(F, M, n, j)
        is_enough, logP = sum_log_probabilities(next_logP, logP)

    if logP == 0: logP *= -1

    return -logP

def log_hyper_probability(F, M, n, j):
    return ( log_c(M, j) + log_c(F - M, n - j) - log_c(F, n) ) / log(10)

def log_c(n, k):
    if (k == n) or not k: return 0
    t = n - k
    if t < k: t = k
    return sum_range(t + 1, n) - sum_factorial(n - t)

def sum_range(minimum, maximum):
    return sum([ log(i) for i in xrange(minimum, maximum + 1) ])

def sum_factorial(n):
    return sum([ log(i) for i in xrange(2, n + 1) ])

def sum_log_probabilities(next_logP, logP):
    if next_logP == 0: return True, logP

    if next_logP > logP:
        common = next_logP
        diff_exponent = logP - common
    else:
        common = logP
        diff_exponent = next_logP - common

    logP = common + ( (log(1 + mpow(10, diff_exponent))) / log(10) )

    # The cumulative summation stops when the increasing is less than 10e-4
    if next_logP - logP > -4: return True, logP

    return False, logP

"""
igraph specific code

"""

def igraph_surprise(g, vc):
    # Parameters F and n are always the same for a certain network
    # calculate F (max. possible number of links in the network)
    F = (g.vcount() * (g.vcount() - 1)) / 2

    # calculate n (actual number of links in the network)
    n = g.ecount()

    # calculate M (max. possible number of intra-community links)
    M = sum([ ((s*(s-1))/2) for s in vc.sizes() ])

    # calculate p (actual number of intra-community links)
    # p is equal to all edges minus community crossing edges
    p = n - vc.crossing().count(True)

    print "F = %s, M = %s, n = %s, p = %s" % (F, M, n, p)
    S = compute_surprise(F, M, n, p)
    return S

def main():
    import igraph
    if len(argv) == 2:
        g = igraph.load(argv[1])
        print g.summary()
        c_im = g.community_infomap()
        s_im = igraph_surprise(g, c_im)
        print c_im.summary()
        print "Infomap: Surprise = %s, Modularity = %s" % (s_im, c_im.q)
        print

        c_fg = g.as_undirected().community_fastgreedy().as_clustering()
        s_fg = igraph_surprise(g, c_fg)
        print c_fg.summary()
        print "Fast Greedy: Surprise = %s, Modularity = %s" % (s_fg, c_fg.q)
        print


        if not g.is_directed():
            c_mls = g.community_multilevel(return_levels=True)
        else:
            c_mls = g.as_undirected().community_multilevel(return_levels=True)
        
        for i, c_ml in enumerate(c_mls):
            s_ml = igraph_surprise(g, c_ml)
            print c_ml.summary()
            print "Multi-Level at Level %s: Surprise = %s, Modularity = %s" % (i, s_ml, c_ml.q)
            print

    elif len(argv) == 3:
        g = igraph.Graph.Read_Edgelist(argv[1])
        
        # load partition data into dict with node id as index
        tmp = dict([ map(int,l.split()) for l in open(argv[2]).readlines() ])

        # re-map to 0 based partition id's
        zmap = dict((v, k) for k,v in enumerate(set(tmp.values())))
        p = dict([(k, zmap[v]) for k,v in tmp.iteritems() ])

        # load partition data into clustering object
        # to re-use igraph functions
        vc = igraph.VertexClustering(g, [ p[v.index] for v in g.vs ])

        print g.summary()
        print vc.summary()

        # get S
        s = igraph_surprise(g, vc)
        print "Surprise = %s, Modularity = %s" % (s, vc.q)
    else:
        print "Usage: %s graph-file" % argv[0]
        print "    or %s network-edge-list partition-file" % argv[0]

if __name__ == '__main__':
    main()

