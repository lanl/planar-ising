
function sigma = kasteleyn_orientation(Gamma,symm_edges,inv_edge)

n = length(Gamma);

if (mod(n,2))
  warning('graph has odd number of nodes, no complete kasteleyn orienation.');
end

m = size(symm_edges,1);

% initially assign random orientations to all edges
sigma = sign(randn(m,1));

% orientation must be anti-symmetric
for k = 1:m
  sigma(inv_edge(k)) = -sigma(k);
end

% construct dual graph
[g,faces,circuits,symm_dual_edges] = dual_graph(Gamma,symm_edges,inv_edge);

f = length(faces); % number of faces

if (f==1)
  return;
end

% take random spanning tree of dual graph

tree_edges = random_spanning_tree(f,symm_dual_edges,inv_edge);

% this returns a directed spanning tree such that:
% 1) edges are oriented to point towards a root node,
% 2) this root node has exactly one incident edge, and 
% 3) tree_edges is listed in leaf-to-root ordering

% adjust orientation on edges cut by dual spanning tree to impose Kasteleyn condition
% on each closed face, proceeding from leaves to root of dual spaning tree...


for k = tree_edges
 
  C = circuits{symm_dual_edges(k,1)};
  s = prod(sigma(C));

  if (s == +1)

    % reverse orientation of cut edge (and its inverse edge)
    sigma(k) = -sigma(k);
    sigma(inv_edge(k)) = -sigma(inv_edge(k));

    % this should fix the face (unless it somehow contains the edge multiple times?)
    % just to be sure...
    s = prod(sigma(C));
    if (s == +1)
      error('something fishy is going on here!');
    end

  end

end

% check if last face is kastelyn (all other faces are)

C = circuits{symm_dual_edges(k,2)};
s = prod(sigma(C));

if (s == +1) 
  % this should happen if (and only if?) odd number of vertices
  warning('one face is not kasteleyn'); 
end
