


edges = [1 2; 1 3; 1 4; 2 3; 3 4; 2 4];
symm_edges = [edges; edges(:,2) edges(:,1)];
inv_edge = [7:12 1:6]';

E = sparse(symm_edges(:,1),symm_edges(:,2),1:12,4,4);


Gamma = cell(4,1);
 
Gamma{1} = [E(1,4) E(1,3) E(1,2)];
Gamma{2} = [E(2,1) E(2,4) E(2,3)];
Gamma{3} = [E(3,2) E(3,1) E(3,4)];
Gamma{4} = [E(4,3) E(4,2) E(4,1)];

subgraph = [1:5 7:11]';


sigma = kasteleyn_orientation(Gamma,symm_edges,inv_edge)


S_1321 = sigma(E(1,3))*sigma(E(3,2))*sigma(E(2,1))
S_1431 = sigma(E(1,4))*sigma(E(4,3))*sigma(E(3,1))
S_12341 = sigma(E(1,2))*sigma(E(2,3))*sigma(E(3,4))*sigma(E(4,1))

[Z_bp,lambda] = planar_cover_infer(zeros(4,1),sigma,symm_edges(:,1),symm_edges(:,2),subgraph)

P_bp = sqrt(1/Z_bp)


