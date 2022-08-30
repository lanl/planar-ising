% Convert exponential parameters of Boltzmann model
% to exponential parameters of equivalent Ising models.
%
% Boltzmann model has variables x = {0,1}
%
% Ising models has variables x = {-1,+1}
%
% Either model can be represented:
% p(x) ~ exp( sum_i theta_i x_i + sum_ij theta_ij x_i x_j )
%
% This code assumes only pairwise interactions are present.

function [theta_nodes,theta_edges,theta_const] = boltz2ising(edges,beta_nodes,beta_edges)

m = length(beta_edges);

theta_edges = 0.25*beta_edges;

theta_nodes = beta_nodes/2.0;

for e = 1:length(beta_edges)
  i = edges(e,1);
  j = edges(e,2);
  theta_nodes(i) = theta_nodes(i) + theta_edges(e);
  theta_nodes(j) = theta_nodes(j) + theta_edges(e);
end
