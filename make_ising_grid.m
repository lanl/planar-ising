% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.


function [mu_nodes,mu_edges,theta_nodes,theta_edges] = make_ising_grid(n)

a = get_graph(n,'grid');
edges = get_edges_from_graph(a);
theta_nodes = rand(n*n,1) - 0.5;
theta_edges = rand(1,size(edges,1)) - 0.5;
[moments_nodes,moments_pairs] = mobius_ising_moments(n*n,edges,theta_nodes,theta_edges)
