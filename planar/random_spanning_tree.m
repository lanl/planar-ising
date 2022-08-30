
function tree_edges = random_spanning_tree(n,symm_edges,inv_edge)

if (n==1)
  tree_edges = [];
  return;
end

m = size(symm_edges,1);

in_tree = zeros(n,1);
parent = zeros(n,1);

root = ceil(n*rand(1));
in_tree(root) = 1;

loc = zeros(m,1);

for k = 1:n-1

  boundary = find(~in_tree(symm_edges(:,1)) & in_tree(symm_edges(:,2)));
  new_edge = boundary(ceil(length(boundary)*rand(1)));

  tree_edges(k) = new_edge;
  new_node = symm_edges(new_edge,1);
  parent(new_node) = new_edge;
  in_tree(new_node) = 1.0;  
  loc(new_edge) = k;

end

% reorient tree so that last node added becomes the root.

% trace path back to root
ij = new_edge;
path = [];
while (1)
  path = [path ij];
  j = symm_edges(ij,2);
  if (j ~= root)
    ij = parent(j);
  else
    break;
  end
end

% reverse orientation of path
inv_path = inv_edge(path)';

tree_edges = [inv_path tree_edges(setdiff([1:n-1],loc(path)))];

% reverse order that tree_edges are listed (to go from leaves to root)
tree_edges = tree_edges(end:-1:1);

