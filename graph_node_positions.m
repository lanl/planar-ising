% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function node_positions = graph_node_positions(n,graph_type)

    if( strcmp(graph_type,'grid') )
        node_positions = zeros(n^2,2);

        for i = 0 : n-1
            for j = 1 : n
                node_positions( n*i + j,: ) = [ j,i+1 ];
            end
        end
    elseif( strcmp(graph_type,'cycle') )
        node_positions = [ 1 : n-1 ]';
        node_positions = [node_positions, ones(n-1,1)];
        node_positions = [ node_positions; 5,5];
    elseif( strcmp(graph_type,'counterexample_modified') )
        node_positions = [ 1 : n-2 ]';
        node_positions = [ ones(n-2,1), node_positions ];
        node_positions = [ 0, 0; node_positions ; 2, 0 ];
    else
        error( 'graph type unknown!' );
    end