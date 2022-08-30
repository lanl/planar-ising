
clear all;

arcs = [1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 5; 3 4; 3 5; 4 5; 5 6];
arcs = [arcs; arcs(:,2) arcs(:,1)];
inv_arc = [12:22 1:11];

E = sparse(arcs(:,1),arcs(:,2),1:20,5,5);

Gamma{1} = [E(1,5) E(1,4) E(1,3) E(1,2)];
Gamma{2} = [E(2,1) E(2,5) E(2,4) E(2,3)];
Gamma{3} = [E(3,2) E(3,1) E(3,5) E(3,4)];
Gamma{4} = [E(4,3) E(4,2) E(4,1) E(4,5)];
Gamma{5} = [E(5,4) E(5,3) E(5,2) E(5,1)];

sigma = kasteleyn_orientation(Gamma,arcs,inv_arc)


subgraph = [1 2 3 4 5 8 10]';
subgraph = [subgraph; inv_arc(subgraph)];


planar_cover_infer(zeros(6,1),sigma,arcs,inv_arc,subgraph)
