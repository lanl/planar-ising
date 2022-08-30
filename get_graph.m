% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function a = get_graph(n,graph_type)

    if( strcmp(graph_type,'cycle') )
        a = sparse(zeros( n ));

        %Constructing a cycle
        for i = 1 : n
            a( i, mod( i , n ) + 1 ) = 1;
            a( mod( i , n ) + 1, i ) = 1;
        end
    elseif( strcmp(graph_type,'grid') )
        a = sparse(zeros(n^2));

        %Constructing a grid
        for i = 0 : n-1
            for j = 1 : n
                if( i ~= n-1 )
                    a( n*i+j, n*i+j+n ) = 1;
                    a( n*i+j+n, n*i+j ) = 1;
                end
                if( j ~= n )
                    a( n*i+j, n*i+j+1 ) = 1;
                    a( n*i+j+1, n*i+j ) = 1;
                end
            end
        end
    elseif( strcmp(graph_type,'complete') )
        a = sparse(ones(n)-eye(n));
    elseif( strcmp(graph_type,'randomouterplanar') )
        a = generate_random_outer_planar_graph(n);
    elseif( strcmp(graph_type,'randomplanar') )
        a = generate_random_planar_graph(n);
    elseif( strcmp(graph_type,'randomgraph') )
        a = generate_random_graph(n,0.2);
    else
        error( 'Graph type unknown!' );
    end
