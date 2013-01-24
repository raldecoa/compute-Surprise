# surprise calculates the parameter Surprise for 4 given parameters
# You should have received this program together with 
# computeSurprise.m (which receives a network and a given partition)
#
# If you use this program, please cite:
#       Aldecoa R, Marín I (2011)
#       Deciphering network community structure by Surprise
#       PLoS ONE 6(9): e24195

# The program receives the four parameters needed for the computation
# of Surprise. It calculates the probability of the partition based on
# a cumulative hypergeometric distribution. All calculations are done
# by using logarithms and other optimizations to avoid underflow problems.
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




surprise <- function(graph, membership) {
  
  # Parameters F and n are always the same for a certain network
  # calculate F (max. possible number of links in the network)
  F <- (vcount(graph)*(vcount(graph)-1))/2
  
  # calculate n (actual number of links in the network)
  n <- ecount(graph)
  
  # calculate M (max. possible number of intra-community links)
  M <- 0
  for(i in unique(membership)){
    s <- sum(membership==i)
    M <- M+((s*(s-1))/2)
  }
  
  # calculate p (actual number of intra-community links)
  p <- 0
  for(i in unique(membership)) {
    # edge ids of community i
    ids <- V(graph)[membership==i]
    p <- p+round(length(E(graph)[ ids %--% ids ]))
  }
  
  S <- surprise_regular(M=M, p=p, F=F, n=n)
  
  return(S)
}

surprise_regular <- function(M, p, F, n) {
  
  min <-  M
  if(n < M) {
    min <-  n
  }
  
  logP <-  logHyperProbability(F=F,M=M,n=n,p=p)
  stop <-  0
  while (stop == 0 && (p < min)) {
    p <-  p+1
    nextLogP <-  logHyperProbability(F,M,n,p)
    result <-  sumLogProbabilities(nextLogP, logP)  
    stop <- result[1]
    logP <- result[2]
  }
  
  return (-logP)
}

logHyperProbability  <-  function(F, M, n, p) {
  
  logC(k=n-p,n1=F-M)
  
  logH <-  logC(k=p,n1=M) + logC(k=n-p,n=F-M) - logC(k=n,n1=F)
  logH <-  logH/log(10)
  return (logH)
}

sumFactorial <- function(n) {
  sum <-  0
  if(n>=2){
    for (i in  2:n) {
      sum <-  sum + log(i)  
    }  
  }
  
  return (sum)
}

sumRange <- function (min, max) {
  sum <-  0
  if(max>=min){
    for (i in  min:max) {
      sum <-  sum + log(i)
    }
  }
  return (sum)
}

logC <- function (k, n1) {
  
  if(k == n1 || !k) {
    x <-  0
  } else {
    t <-  n1 - k
    if(t < k) {
      t <-  k
    }
    
    x <-  sumRange(min=t+1, max=n1) - sumFactorial(n1 - t)
  }
  
  return (x)
}



sumLogProbabilities <- function(nextLogP, logP) {
  stop <- 0
  
  if(nextLogP == 0) {
    stop <- 1
  } else {
    stop <- 0
    if(nextLogP > logP) {
      common <-  nextLogP
      diffExponent <-  logP - common
    }else {
      common <-  logP
      diffExponent <-  nextLogP - common
    }
    
    logP <-  common + ((log(1 + 10^diffExponent)) / log(10))
    
    if(nextLogP - logP > -4) {
      stop <- 1
    }
  }
  
  return (c(stop, logP))
}