
% Convert exponential parameters of Ising model
% to exponential parameters of equivalent Boltzmann model.
%
% Boltzmann model has variables x = {0,1}
% Ising models has variables x = {-1,+1}
%
% Either model can be represented:
% p(x) ~ exp( sum_i theta_i x_i + sum_ij theta_ij x_i x_j )
%
%
% This code assumes only pairwise interactions are present.

function [beta_nodes,beta_edges,beta_const] = ising2boltz(edges,theta_nodes,theta_edges)

m = length(theta_edges);

beta_edges = 4.0*theta_edges;

beta_nodes = theta_nodes;
for e = 1:m
  i = edges(e,1);
  j = edges(e,2);
  beta_nodes(i) = beta_nodes(i) - theta_edges(e);
  beta_nodes(j) = beta_nodes(j) - theta_edges(e);
end
beta_nodes = 2.0*beta_nodes;

beta_const = sum(theta_edges) - sum(theta_nodes);
