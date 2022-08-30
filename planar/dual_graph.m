
% function [g,faces,circuits,symm_dual_edges,inv_dual_edge] = dual_graph(Gamma,symm_edges,inv_edge)
%
% SYNOPSIS: construct dual graph of an embedded surface graph.
%
% INPUT:
% Gamma: specifies embedded graph by ccw-cyclic ordering of each node's incident edges
% symm_edges: the symmetric edge set of graph
% 
% OUTPUT:
% g: the genus of  embedded graph with circular ordering gamma.
% faces: cell array of the faces of the graph (cw ordering of vertices)
% symm_dual_edges: the symmetric edge set of dual graph (linking faces of graph) 
%

function [g,faces,circuits,symm_dual_edges] = dual_graph(Gamma,symm_edges,inv_edge)

n = length(Gamma);
m = size(symm_edges,1);

symm_dual_edges = zeros(m,2);
faces = {};
f = 1; % face counter

for k = 1:m

  if (symm_dual_edges(k,2) == 0)

    % walk around the face to left of edge k
    A = [];
    C = [];

    ij = k;

    while (1)

      i = symm_edges(ij,1);
      j = symm_edges(ij,2);
      ji = inv_edge(ij);

      %fprintf('%d -> %d\n',i,j);
     
      symm_dual_edges(ij,2) = f;
      symm_dual_edges(ji,1) = f;
      
      A = [A i];
      C = [C ij];

      % move to next edge
      gamma_j = Gamma{j};
      d_j = length(gamma_j);
      t = find(gamma_j==ji);
      if (t>1) 
        t=t-1;
      else
        t=d_j;
      end
      ij = gamma_j(t);

      % stop when the circuit is closed
      if (ij==k) 
        break; 
      end 

    end

    faces{f} = A;
    circuits{f} = C;

    f = f+1;

  end

end

f = f-1; % number of faces

g = 1+(m/2-n-f)/2; % compute genus using euler's formula
