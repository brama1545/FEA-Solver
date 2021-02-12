clc;
clear;

nodes = 4; %INPUT # of nodes
edges = 3; %INPUT # of elements in sys

nodeFile = "Spring.txt"; %See doc for formatting
BCFile = "boundaries.txt"; %See doc for formatting

nodeFileNo = fopen(nodeFile, 'r');
BCFileNo = fopen(BCFile, 'r');
input = fscanf(nodeFileNo, "%f", [3, edges]);
boundaries = fscanf(BCFileNo, "%f", [3, Inf]);
fclose(nodeFileNo);
fclose(BCFileNo);

nodeList = struct('ID', {1:nodes}, 'Edges', {1:nodes*[]}, 'k', [], 'BC', [], 'Load', [], 'Theta', []);
%nodeList = repmat(nodeStruct, 1, nodes); %List of node structs with props

STIFFNESS = zeros(nodes); %Global stiffness matrix

for edge = 1:edges
    node1 = input(1, edge);
    node2 = input(2, edge);
    
    nodeList(node1).Edges(end + 1) = node2;
    nodeList(node2).Edges(end + 1) = node1;
    
    nodeList(node1).k(end + 1) = input(3, edge);
    nodeList(node2).k(end + 1) = input(3, edge);
    
    STIFFNESS(node1, node2) = STIFFNESS(node1, node2) - input(3, edge);
    STIFFNESS(node2, node1) = STIFFNESS(node2, node1) - input(3, edge);
    STIFFNESS(node1, node1) = STIFFNESS(node1, node1) + input(3, edge);
    STIFFNESS(node2, node2) = STIFFNESS(node2, node2) + input(3, edge);
end

BCS = zeros(nodes, 1);
SStiffness = STIFFNESS;
fixedNodes = [];
for bound = 1:length(boundaries)
    if ~boundaries(1, bound)
        fixedNodes = [fixedNodes, boundaries(2, bound)];
        BCS(boundaries(2, bound)) = boundaries(3, bound);
    end
end

annihilator = STIFFNESS * BCS;

inForce = zeros(nodes, 1);
for force = 1:length(boundaries)
    if boundaries(1, bound)
        inForce(boundaries(2, bound)) = boundaries(3, bound);
    end
end

inForce = inForce - annihilator;
for node = fixedNodes
    inForce(node) = [];
    SStiffness(boundaries(2, bound), :) = [];
    SStiffness(: ,boundaries(2, bound)) = [];
end
displacement = inForce\SStiffness;


