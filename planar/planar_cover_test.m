

% build 2D grid of dimensions w by w

w = 4;
n = w*w;
m = 2*w*(w-1)+2*(w-1)^2;

grid = reshape([1:n],w,w);

edges = zeros(m,2);

sub1 = zeros(m,1);
sub2 = zeros(m,1);

k=1;

% vertical edges
for i = 1:w-1
  for j = 1:w

    edges(k,1) = grid(i,j);
    edges(k,2) = grid(i+1,j);    
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
    sub1(k)=mod(k,2);
    sub2(k)=mod(k+1,2);
    k = k+1;
  end
end

for i = 1:w-1
  for j = 1:w-1
    edges(k,1) = grid(i+1,j);
    edges(k,2) = grid(i,j+1);    
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

% Build information matrix H

r = .1 * ones(2*m,1);
H_arcs = -r;
H_diag = ones(n,1);


H = diag(H_diag) + sparse(arcs(:,1),arcs(:,2),H_arcs,n,n);
log_Z = -log(det(H)) 

H1 = diag(H_diag) + sparse(arcs(arcs1,1),arcs(arcs1,2),H_arcs(arcs1),n,n);
log_Z1 = -log(det(H1))

H2 = diag(H_diag) + sparse(arcs(arcs2,1),arcs(arcs2,2),H_arcs(arcs2),n,n);
log_Z2 = -log(det(H2))

H12 = diag(H_diag) + sparse(arcs(arcs12,1),arcs(arcs12,2),H_arcs(arcs12),n,n);
log_Z12 = -log(det(H12))


log_Z_combine = log_Z1 + log_Z2 - log_Z12


% test planar cover method

[Z_pc1,lambda] = planar_cover_infer(H_diag,H_arcs,arcs(:,1),arcs(:,2),arcs1);
log_Z_pc1 = log(Z_pc1)

[Z_pc2,lambda] = planar_cover_infer(H_diag,H_arcs,arcs(:,1),arcs(:,2),arcs2);
log_Z_pc2 = log(Z_pc2)

[Z_pc12,lambda] = planar_cover_infer(H_diag,H_arcs,arcs(:,1),arcs(:,2),arcs12);
log_Z_pc12 = log(Z_pc12)

log_Z_pc_combine = log_Z_pc1 + log_Z_pc2 - log_Z_pc12


lambda




