% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.


% Jason K. Johnson
% LANL, Sep 2008.

% OVERVIEW: Kac-Ward determinant calculation for partition sum of
% planar zero-field Ising model.

function [Z,M,H] = kacward(points,edges,K,require_hessian)

% INPUTS
% points: vertex locations in the plane
% edges: matrix specifying pairs of vertices
% K: edge weights of Ising model
% If the last edge has 0 weight, moments will be calculated for that
% If any other edge has 0 weight, the moment will be NaN

% OUTPUT: Z, the partition sum of Ising model.
%         M, moments on the edges of the Ising model.

n = size(points,1);
m = size(edges,1);


% build list of arcs (split each undirected edge into 2 directed edges)

redges = edges;
redges(:,1) = edges(:,2);
redges(:,2) = edges(:,1);
arcs = [edges; redges];

% build lists of edges into and out of each vertex

arcs_in = cell(n,1);
arcs_out = cell(n,1);

for v = 1:n
  arcs_in{v} = find(arcs(:,2)==v)';
  arcs_out{v} = find(arcs(:,1)==v)';
end

% build arc-to-arc transition matrix

A = zeros(2*m,2*m);
w = tanh(K);
w = [w, w];
w_matrix = diag(w);

for v = 1:n
  for arc1 = arcs_in{v}
    for arc2 = arcs_out{v}

      if (arcs(arc1,2)~=v | arcs(arc2,1)~=v)
        abort;
      end

      s = arcs(arc1,1);
      t = arcs(arc2,2);
 
      if (s~=t) % non-backtracking walks

        p1 = points(s,:)';
        p2 = points(v,:)';
        p3 = points(t,:)';

        % compute phase factor based on turning angle
        r1 = p2-p1;
        r2 = p3-p2;
        phi1 = atan2(r1(2),r1(1));
        phi2 = atan2(r2(2),r2(1));
        phi = mod(phi2-phi1+pi,2*pi)-pi;

        A(arc2,arc1) = exp(sqrt(-1)*phi/2);
      end

    end
  end
end


A = sparse(A);
I = speye(2*m);
AW = sparse(A*w_matrix);

% compute partition sum using determinant calculation

zeta = det(I-AW);

if (abs(imag(zeta)/real(zeta))>1e-6)
  error('determinant not real.');
end

if (real(zeta)<0.0)
  error('determinant negative.');
end

Z = 2^n * prod(cosh(K)) * sqrt(real(zeta));

%F = log(2) + (sum(log(cosh(K))) + 0.5 * log(real(zeta)))/n;
%F = (sum(log(cosh(K)))+0.5*log(real(zeta)))/n;

B = inv(I-AW);
U = B*A;

M = zeros(1,m);
H = zeros(m);

if( require_hessian )
    for index = 1 : m
        M(index) = real( w(index) - 0.5* ( (1-w(index)^2) * ( U(index,index)+U(m+index,m+index) ) ) );
        for index_inner = 1 : m
            H(index,index_inner) = -0.5 * real( (1-w(index)^2) * (U(index,index_inner)*U(index_inner,index)+U(index,m+index_inner)*U(m+index_inner,index)+U(m+index,index_inner)*U(index_inner,m+index)+U(m+index,m+index_inner)*U(m+index_inner,m+index)) * (1-w(index_inner)^2) );
            if( index == index_inner )
                %H(index,index_inner) = H(index,index_inner) + (1-w(index)^2);
                H(index,index) = 1 - M(index)^2;
            end
        end
    end
else
    for index = 1 : m
        M(index) = real( w(index) - 0.5 * ( (1-w(index)^2) * (U(index,index)+U(m+index,m+index)) ) );
    end
end