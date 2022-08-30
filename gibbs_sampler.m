% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

% samples = gibbs_sampler(g, nsamples, beta=1)
%
% g is a weighted undirected graph.
% samples is a nsamples x N matrix of +1/-1 values
% optional argument, beta

function samples = gibbs_sampler(g, nsamples, beta)

if (nargin <= 3)
    beta = 1;
end

[N, ~] = size(g);

% Initialize
samples = zeros(nsamples, N);
currState = randi(2, [1, N]);
currState(currState==2) = -1;

% Draw samples
for t = 1:nsamples
    % randomize order of variables
    nodes = randperm(N);
    for i = 1:N
        node = nodes(i);
        probPos = (1 + exp(-2 * beta * sum(g(node,:) .* currState))) ^ (-1);
        sigma = rand > probPos;
        if (~sigma)
            currState(node) = -1;
        else
            currState(node) = 1;
        end
    end
    samples(t, :) = currState;
end
