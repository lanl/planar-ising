
% build 2D grid of dimensions w by w

w = 4;
n = w*w;
m = 2*w*(w-1)+2*(w-1)^2;

grid = reshape([1:n],w,w);

edges = zeros(m,2);
sigma = zeros(m,1);


sub1 = zeros(m,1);
sub2 = zeros(m,1);

k=1;

% vertical edges
for i = 1:w-1
  for j = 1:w

    edges(k,1) = grid(i,j);
    edges(k,2) = grid(i+1,j);    
    sigma(k) = +1;

    sub1(k)=1;
    sub2(k)=1;

    k = k+1;

  end
end


% horizontal edges
for i = 1:w
  for j = 1:w-1

    edges(k,1) = grid(i,j);
    edges(k,2) = grid(i,j+1);    
    sigma(k) = 2*mod(i,2)-1;

    sub1(k)=1;
    sub2(k)=1;
    k = k+1;

  end
end

% "crossed" diagonal edges (non-planar)
for i = 1:w-1
  for j = 1:w-1

    edges(k,1) = grid(i,j);
    edges(k,2) = grid(i+1,j+1);    
    sigma(k) = 2*mod(i+j,2)-1; 

    sub1(k)=mod(k,2);
    sub2(k)=mod(k+1,2);
    k = k+1;

  end
end

for i = 1:w-1
  for j = 1:w-1

    edges(k,1) = grid(i+1,j);
    edges(k,2) = grid(i,j+1);    
    sigma(k) = 2*mod(i-j,2)-1; 

    sub1(k)=mod(k,2);
    sub2(k)=mod(k+1,2);
    k = k+1;

  end
end

arcs = [edges; edges(:,2) edges(:,1)];

sub1 = [sub1; sub1];
sub2 = [sub2; sub2];
sub12 = sub1 .* sub2;

arcs1 = find(sub1);
arcs2 = find(sub2);
arcs12 = find(sub12);

sigma = [sigma -sigma];


% build oriented (skew-symmetric) adjacency matrix of the graph

A_diag = zeros(n,1);
A_arcs = sigma;


A = sparse(arcs(:,1),arcs(:,2),sigma,n,n);

A1 = sparse(arcs(arcs1,1),arcs(arcs2,2),sigma(arcs1),n,n);
P1 = sqrt(det(A1))

A2 = sparse(arcs(arcs2,1),arcs(arcs2,2),sigma(arcs2),n,n);
P2 = sqrt(det(A2))

A12 = sparse(arcs(arcs12,1),arcs(arcs12,2),sigma(arcs12),n,n);
P12 = sqrt(det(A12))

P_add = P1 + P2 - P12

P_mul = P1 * P2 / P12


[Z1,lambda] = planar_cover_infer(A_diag,A_arcs,arcs(:,1),arcs(:,2),arcs1)

P1_pc = sqrt(1/Z1)

[Z2,lambda] = planar_cover_infer(A_diag,A_arcs,arcs(:,1),arcs(:,2),arcs2)

P2_pc = sqrt(1/Z2)

[Z12,lambda] = planar_cover_infer(A_diag,A_arcs,arcs(:,1),arcs(:,2),arcs12)

P12_pc = sqrt(1/Z12)

P_pc_add = P1_pc + P2_pc - P12_pc

P_pc_mul = P1_pc * P2_pc / P12_pc


