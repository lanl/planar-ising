% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function [ Z, moments, H ] = kacward_allmoments( a, edges, theta_vec )

    n = length(a);
    node_positions = chrobak_payne_straight_line_drawing( sparse(a) );

    [ Z, moments_vec, H ] = kacward( node_positions, edges, theta_vec, true );

    moments = zeros( n );
    
    for index = 1 : size(edges,1)
        moments(edges(index,1),edges(index,2)) = moments_vec(index);
    end
    
    [ ci, sizes ] = components(a);
    
    nontrivial_comp = find(sizes~=1);
    
    for index = 1 : length(nontrivial_comp)
        compnum = nontrivial_comp(index);
        nodes_in_comp = find(ci == compnum);
        num_nodes_in_comp = sizes(compnum);
        a_subcomp = a(nodes_in_comp,nodes_in_comp);
        theta_vec_subcomp = zeros(1,sum(sum(a_subcomp))/2);
        edges_subcomp = zeros(sum(sum(a_subcomp))/2,2);
        index_count = 1;
        for index_inner = 1 : length(theta_vec)
            if( ci(edges(index_inner,1)) == compnum && ci(edges(index_inner,2)) == compnum )
                i_subcompindex = find(nodes_in_comp == edges(index_inner,1));
                j_subcompindex = find(nodes_in_comp == edges(index_inner,2));
                edges_subcomp(index_count,:) = [i_subcompindex,j_subcompindex];
                theta_vec_subcomp(index_count) = theta_vec(index_inner);
                index_count = index_count + 1;
            end
        end
        if (size(edges_subcomp,1) < num_nodes_in_comp)
            for i = 1 : num_nodes_in_comp
                for j = i+1 : num_nodes_in_comp
                    if( a_subcomp(i,j) ~= 0 )
                        continue;
                    end
                    i_actual = nodes_in_comp(i);
                    j_actual = nodes_in_comp(j);
                    a_subcomp(i,j) = 1;
                    a_subcomp(j,i) = 1;
                    if( ~test_planar_graph( a_subcomp ) )
                        moments(i_actual,j_actual) = NaN( 'double' );
                    else
                        node_positions = chrobak_payne_straight_line_drawing( a_subcomp );
                        [ ~, moments_vec, ~ ] = kacward( node_positions, [ edges_subcomp; i, j ], [ theta_vec_subcomp, 0 ], false );
                        moments(i_actual,j_actual) = moments_vec( length(moments_vec) );
                    end
                    a_subcomp(i,j) = 0;
                    a_subcomp(j,i) = 0;
                end
            end
        else
            complete_graph = sparse(ones(num_nodes_in_comp,num_nodes_in_comp) - eye(num_nodes_in_comp));
            a_complement = complete_graph - a_subcomp;
            edges_not_present = get_edges_from_graph(a_complement);
            edges_not_present_planar = [];

            for index = 1 : size(edges_not_present,1)
        %         disp(index);
                i = edges_not_present(index,1);
                j = edges_not_present(index,2);
                i_actual = nodes_in_comp(i);
                j_actual = nodes_in_comp(j);
                a_subcomp(i,j) = 1;
                a_subcomp(j,i) = 1;
                if( ~test_planar_graph( a_subcomp ) )
                    moments(i_actual,j_actual) = NaN( 'double' );
                else
                    edges_not_present_planar = [ edges_not_present_planar ; edges_not_present(index,:) ];
                end
                a_subcomp(i,j) = 0;
                a_subcomp(j,i) = 0;
            end

            while size(edges_not_present_planar,1) > 0
                a_new = a_subcomp;
                edges_new = edges_subcomp;
                theta_vec_new = theta_vec_subcomp;
                index = 1;
                while index <= size(edges_not_present_planar,1)
                    i = edges_not_present_planar(index,1);
                    j = edges_not_present_planar(index,2);
                    a_new(i,j) = 1;
                    a_new(j,i) = 1;
                    if( ~test_planar_graph( a_new ) )
                        a_new(i,j) = 0;
                        a_new(j,i) = 0;
                        index = index+1;
                    else
                        edges_new = [ edges_new ; i,j ];
                        theta_vec_new = [ theta_vec_new , 0 ];
                        edges_not_present_planar = [ edges_not_present_planar(1:index-1,:);edges_not_present_planar(index+1:size(edges_not_present_planar,1),:) ];
                    end
                end
                node_positions = chrobak_payne_straight_line_drawing( a_new );
                [ ~, moments_vec, ~ ] = kacward( node_positions, edges_new, theta_vec_new, false );

                for index = size(edges_subcomp,1)+1 : size(edges_new,1)
                    i = edges_new(index,1);
                    j = edges_new(index,2);
                    i_actual = nodes_in_comp(i);
                    j_actual = nodes_in_comp(j);
                    moments(i_actual,j_actual) = moments_vec( index );
                end
            end
        end
    end
    moments = moments + moments';