computeSurprise <-  function(networkFile, partitionFile) {

# computeSurprise calculates the Surprise of the partition of a network
# This program needs surprise.m to work properly. You should have 
# received surprise.R along with this program.
# 
# If you use this program, please cite:
#       Aldecoa R, Marín I (2011)
#       Deciphering network community structure by Surprise
#       PLoS ONE 6(9): e24195
# 
# The program receives two input files:
#  - networkFile: A network represented by a list of links (pairs of nodes)
#                 Each line contains two nodes separated by a 'tab'
#  - partitionFile: Describes a given partition of the network
#                 Each line contains a node and the community to which it
#                 is assigned, separated by a 'tab'.
#                 (The partition identifier must be a number)
#  
# ** Two toy examples of these files are included within this folder 
#    (network.pairs and partition.part)
#
# The input graph is represented based on the graph class in the igraph package
# http://cran.r-project.org/web/packages/igraph
# 
# Copyright (C) 2012 Rodrigo Aldecoa and Ignacio Marín
# 
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
# 
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
#   Contact info: Rodrigo Aldecoa <raldecoa@ibv.csic.es>
  

   library(igraph)
    
   # read network and partition
   network <- read.table(networkFile)
   part <- read.table(partitionFile)
   g <- graph.edgelist(as.matrix(network), directed=FALSE)

   part <- part[order(part$V1),] # sort table by node index
   mem <- part[,2] # get communities from table
   surprise(graph=g, membership=mem)  
}


computeSurprise(networkFile="network.pairs",
                partitionFile="partition.part")
