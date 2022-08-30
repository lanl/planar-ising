% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function logz = mrf_logz(nodePot, edgePot, edgeEnds)
% LL = MRF_LOGLIKELIHOOD(Y, NODEPOT, EDGEPOT, EDGEENDS)
%
% Calculates the log partition function of a MRF.
% logz           = log partition function of given model.
% nodePot(j,v)   = potential function for node j with value v.
% edgePot(u,v,e) = potential function for edge e with node1 value u and
%                  node2 value v, maxvalue of e is number of non-zero
%                  edges.
% edgeEnds(e,:)  = pair of nodes on ends of edge e.

nNodes = size(nodePot, 1);
nEdges = size(edgeEnds, 1);

% Convert to log space
nodeF = log(nodePot);
edgeF = log(edgePot);
ll = 0;


energy1 = 0;
energy2 = 0;
entropy1 = 0;
entropy2 = 0;
