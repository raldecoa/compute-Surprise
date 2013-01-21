function S = computeSurprise(networkFile, partitionFile)

% computeSurprise calculates the Surprise of the partition of a network
% This program needs surprise.m to work properly. You should have 
% received surprise.m along with this program.
%
% If you use this program, please cite:
%       Aldecoa R, Marín I (2011)
%       Deciphering network community structure by Surprise
%       PLoS ONE 6(9): e24195

% The program receives two input files:
%  - networkFile: A network represented by a list of links (pairs of nodes)
%                 Each line contains two nodes separated by a 'tab'
%  - partitionFile: Describes a given partition of the network
%                 Each line contains a node and the community to which it
%                 is assigned, separated by a 'tab'.
%                 (The partition identifier must be a number)
%  
% ** Two toy examples of these files are included within this folder 
%    (network.pairs and partition.part)

% Copyright (C) 2012 Rodrigo Aldecoa and Ignacio Marín
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Contact info: Rodrigo Aldecoa <raldecoa@ibv.csic.es>



% Surprise parameters are F, M, n and p

% READ PARTITION FILE
pFile = fopen(partitionFile, 'rt');
tmp = textscan(pFile, '%s');
fclose(pFile);
nobjects = size(tmp{1},1);
nnodes = nobjects/2;
nodeLabels = containers.Map();
commLabels = [];
for i = 1:2:nobjects
    nodeLabels(tmp{1}{i}) = idivide(i, int32(2), 'floor')+1;
    commLabels = [commLabels str2num(tmp{1}{i+1})];
end

ncomm = numel(unique(commLabels'));
commsMap = containers.Map(unique(commLabels), (1:ncomm));

% Create partition structures
comm = cell(1,ncomm);
node2partition = zeros(1,nnodes);
for i=1:nnodes
    c = commsMap(commLabels(i));
    comm{c} = [comm{c} i]; % Add node to community
    node2partition(i) = c; % Assign community to a node
end

% READ NETWORK FILE
nwFile = fopen(networkFile, 'rt');
links = textscan(nwFile, '%s');
nlinks = size(links{1},1)/2;

n = nlinks; % Parameter n (number of links in the network)

% Parameter F (maximum possible number of links in the network)
F = nnodes * (nnodes-1) / 2;

%Parameter M (maximum possible number of intra-community links)
M=0;
for i = 1:ncomm
    commSize = size(comm{i},2);
    M = M + (commSize * (commSize-1) / 2);
end


% Create the adjacency list of the graph
% and obtain parameter p (number of intra-community links
adjacencyList = cell(1,nnodes);
p = 0;
for i = 1:2:(nlinks*2)
    nodeA = nodeLabels(links{1}{i});
    nodeB = nodeLabels(links{1}{i+1});
    if(nodeA > nodeB)
        tmp = nodeB;
        nodeB = nodeA;
        nodeA = tmp;
        clear tmp;
    end
    adjacencyList{nodeA} = [adjacencyList{nodeA} nodeB];
    if(node2partition(nodeA) == node2partition(nodeB))
        p = p + 1;
    end
end

% COMPUTE SURPRISE
S = surprise(F, M, n, p);