% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function test_kacward()

%Write the graph here

n = 3;
a = get_graph(n,'grid');
n = size(a,1);
edges = get_edges_from_graph(a);
theta = randn(1,size(edges,1));

[Z_brute,M_brute,H_brute] = partition_moments_hessian(n,edges,theta);

[Z,M,H] = kacward_allmoments( a, edges, theta );

Z-Z_brute
M-M_brute
H-H_brute
