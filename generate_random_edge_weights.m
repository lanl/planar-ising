% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

% g = generate_random_edge_weights(a, threshold)
%
% a is an undirected binary graph
% g is an undirected weights graph of the same size
% Weights are chosen uniformly random on the interval [-1,1], but the
% absolute value must be greater than given threshold (default 0.05).

function g = generate_random_edge_weights(a, threshold)

if (nargin <= 2)
    threshold = 0.05;
end

[n, ~] = size(a);

% need weights
tmp = -0.95 + 1.9 * rand(n);
tmp = tmp + sign(tmp) * threshold;
tmp = triu(tmp,1);
tmp = tmp + tmp';
g = a .* tmp;
clear tmp
