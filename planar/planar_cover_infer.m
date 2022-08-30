
function [Z,lambda] = planar_cover_infer(H_diag,H_arcs,arcs,inv_arc,subgraph)

arc_src = arcs(:,1);
arc_tgt = arcs(:,2);

n = length(H_diag);  % # nodes in graph
m = length(H_arcs);  % # arcs in graph (twice # of undirected edges)
m_sub = length(subgraph); % # of arcs in subgraph

% build (complete) H matrix
H = diag(H_diag) + sparse(arc_src,arc_tgt,H_arcs,n,n);

% identify "cut" edges
cut_arcs = setdiff([1:m],subgraph);
m_cut = length(cut_arcs);

arc2cut = zeros(m,1);
arc2cut(cut_arcs) = [1:m_cut];
inv_cut_arc = arc2cut(inv_arc(cut_arcs));

% build H matrix of subgraph
H_0 = diag(H_diag) + sparse(arc_src(subgraph),arc_tgt(subgraph),H_arcs(subgraph),n,n);

% initialize vector of messages on cut arcs
lambda = ones(m_cut,1);

% run the iterative message-passing algorithm

for t = 1:500

    H_sub = H_0;

    % add incoming messages to subgraph
    for k = 1:m_cut
      ij = cut_arcs(k);
      j = arc_tgt(ij);
      H_sub(j,j) = H_sub(j,j) + lambda(k);
    end

    % inference within subgraph
    K_sub = inv(H_sub);

    % update messages

    lambda_old = lambda;

    for k = 1:m_cut
      ij = cut_arcs(k);
      i = arc_src(ij);
      j = arc_tgt(ij);
      lambda_ji = lambda_old(inv_cut_arc(k));
      % update lambda_ij
      lambda(k) = - H(j,i) * (1/K_sub(i,i)-lambda_ji)^(-1) * H(i,j);
    end

    delta = max(abs(lambda-lambda_old));

    %fprintf('t=%d, delta=%g\n',t,delta);

    if (delta<1e-8) 
      break;
    end

end

% compute planar-cover estimate of Z = det(K), where K = inv(H)

H_sub = H_0;

% add incoming messages to subgraph

for k = 1:m_cut
  ij = cut_arcs(k);
  j = arc_tgt(ij);
  H_sub(j,j) = H_sub(j,j) + lambda(k);
end

% inference in subgraph
K_sub = inv(H_sub);
Z = det(K_sub);

% correction for cut edges

for k = 1:m_cut

  ij = cut_arcs(k);
  i = arc_src(ij);
  j = arc_tgt(ij);

  if (i<j)
    K_i = K_sub(i,i);
    K_j = K_sub(j,j);
    lambda_ij = lambda(k);
    lambda_ji = lambda(inv_cut_arc(k));
    H_ij = H_arcs(ij); 
    H_ji = H_arcs(inv_arc(ij));
    H_cut = [(1/K_i-lambda_ji)  H_ij;  H_ji (1/K_j-lambda_ij)];
    m_ij = 1.0 / (det(H_cut)*K_i*K_j);
    Z = Z * m_ij;
  end

end

