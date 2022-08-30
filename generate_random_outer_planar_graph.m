% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function a = generate_random_outer_planar_graph(n)

    a = sparse(zeros(n+1));
    
    a(1,:) = a(1,:)+1;
    a(:,1) = a(:,1)+1;
    a(1,1) = 0;
    
    count = n;
    iter  = n;
    
    for i = 1 : n
        new_edges = randi([2,n+1],count,2);
        for j = 1 : size(new_edges,1)
            node1 = new_edges(j,1);
            node2 = new_edges(j,2);
            if( node1 ~= node2 && a(node1,node2) == 0 )
                a(node1,node2) = 1;
                a(node2,node1) = 1;
                if( ~test_planar_graph(a) )
                    a(node1,node2) = 0;
                    a(node2,node1) = 0;
                end
            end
        end
    end
    
    a = a( 2:n+1,2:n+1);