% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.


function x = graph_embed(edges,x0,it,step)

  n = max(max(edges));
  m = size(x0,1);

  E = zeros(n,n);

  nbrs = cell(n,1);
  for k = 1:m
    i = edges(k,1);
    j = edges(k,2);
    E(i,j) = k;
    E(j,i) = k;
    nbrs{i} = union(nbrs{i},j);
    nbrs{j} = union(nbrs{j},i);
  end

  x = x0;
  g = zeros(size(x));
  lambda = .1;

  for t = 1:it

    g(:) = 0;

    % calculate gradient
    for i = 1:n
      for j = nbrs{i}
        delta = x(j,:) - x(i,:);
	g(i,:) = g(i,:) + (norm(delta)-1) * delta;
      end
    end

    % update
    x = x + step*g;

  end

return;

