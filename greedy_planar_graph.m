% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function [undirected_graph,algo_progress_profile,edge_addition_order,edges_profile,theta_profile] = ...
    greedy_planar_graph( n, moments, aux, numinit, meanweight, imppeople, undirectedgraph_initial, ...
    pruneonthego, addpruneonthego )

    % Default arguments
    if (nargin < 10)
        addpruneonthego = false;
        if (nargin < 9)
            pruneonthego = false;
            if (nargin < 8)
                undirectedgraph_initial = zeros(n);
                if (nargin < 7)
                    imppeople = false;
                    if (nargin < 6)
                        meanweight = 1;
                        if (nargin < 5)
                            numinit = 0;
                            if (nargin < 4)
                                aux = false;
                            end
                        end
                    end
                end
            end
        end
    end

    %improvement_threshold = 3e-4;%1e-5;%1e-10;%3e-4;;
    improvement_threshold = 0;
    lowest_improvement_threshold = inf; %begin with 0
    edge_importance_scaledownfactor = 0.5;
    edge_importance_scaledownfactor_constant = 0.9;
	mutual_info = zeros( n );

	%We assume all the random variables are zero mean and construct a lower triangular matrix
	for i = 1 : n
		for j = i + 1 : n
			mutual_info( i, j ) = 0.5 * ( ( ( 1 + moments( i, j ) ) * log( 1 + moments( i, j ) ) ) + ( ( 1 - moments( i, j ) ) * log( 1 - moments( i, j ) ) ) );
			mutual_info( j, i ) = mutual_info( i, j );
		end
	end

% 	% Run Kruskal's algorithm on the negative mutual information matrix and get the (lower triangular) spanning tree as output
% 	[ graph ] = kruskal_mst( sparse(-mutual_info) );

%   Adding the first edge

    [y_max, j_max] = max(mutual_info);
    [~, i_max] = max(y_max);
% 	graph = ( graph ~= 0 );
%     graph = sparse(zeros(n));
    
    % warned that sparse graph is slow if adding non-zero entries later
    %graph = sparse(undirectedgraph_initial);
    graph = undirectedgraph_initial;
    
%     %adding all edges with node 1 for self potential graphical models
%     graph(:,1)=1;
%     graph(1,:)=1;
%     graph(1,1)=0;

    n_iter = 3 * n - 6;
    edges_info = zeros(n_iter, 3);
    nll_info = zeros(n_iter + 1, 1);
    nllimprovement_lowerbound_info = zeros(n_iter, 1);

    if( ~imppeople && ( ~aux || numinit == 0 ) && sum(sum(graph)) == 0 )
        %adding the edge with the maximum mutual information
        graph(i_max, j_max(i_max)) = 1;
        graph(j_max(i_max), i_max) = 1;
        edges_info(1,:) = [i_max, j_max(i_max), 1];
        nllimprovement_lowerbound_info(1) = mutual_info(i_max,j_max(i_max));
    elseif( imppeople  && sum(sum(graph)) == 0 )
        graph(1,3) = 1;
        graph(3,1) = 1;
        graph(1,41) = 1;
        graph(41,1) = 1;
        graph(1,51) = 1;
        graph(51,1) = 1;
    elseif( sum(sum(graph)) == 0 )
        means = moments(1,2:n);
        [~, indicesorder] = sort(means,'descend');
        for count = 1 : numinit
            graph(1,indicesorder(count+1)) = 1;
            graph(indicesorder(count+1),1) = 1;
        end
    end
    undirected_graph = graph;
    graph = tril( graph );
    complete_graph = ones( n ) - eye( n );
	remaining_edges = tril( complete_graph ) - graph;
    edge_addition_order = zeros(n_iter, 2);
    edge_addition_order(1, :) = [i_max,j_max(i_max)];
    theta_profile = cell(n_iter + 1,1);
    edges_profile = cell(n_iter + 1,1);
    %edges_already_added = [ get_edges_from_graph(undirected_graph), -ones(sum(sum(graph)),1) ];
%     algo_progress_profile{1} = edges_already_added;

%     theta_profile{1} = [];
%     edges_profile{1} = [];
    iter = 1;
%     continueiterations = true;
%     while( continueiterations )

        while( sum( sum( remaining_edges ) ) > 0 )
            iter = iter+1;
    %         if iter == 16
    %             iter = 16;
    %         end
            tmp = ['    Iter ' int2str(iter) ' at time ' datestr(clock)];
            %disp('iter');
            %disp(iter);
            disp(tmp)
            %save 12outerplanar_debugged_traintest_intermediateresults.mat undirected_graph -MAT;
            sparse_graph = sparse(undirected_graph);
            [moments_current_graph,edges_current_graph,theta_current_graph] = get_projection( sparse_graph, moments );

            node_positions = chrobak_payne_straight_line_drawing(sparse_graph);
            [ partition_function, ~ ] = kacward( node_positions, edges_current_graph, theta_current_graph, false );
            loglikelihood = -log(partition_function);
            for edgenumber = 1 : size(edges_current_graph,1)
                loglikelihood = loglikelihood + ( theta_current_graph(edgenumber) * ...
                    moments(edges_current_graph(edgenumber,1),edges_current_graph(edgenumber,2)) );
            end
            nll_info(iter) = loglikelihood;

            if(addpruneonthego)
                node_positions = chrobak_payne_straight_line_drawing( undirected_graph );
                [~,~,~,moments_edges_removed] = kacward_momentsedgeremoved(node_positions,edges_current_graph,theta_current_graph,false);

                maxaicimprovement = -1;
                maxaicimprovementedge = [ -1 -1];
                for inneriter = 1 : length(theta_current_graph)
                    m_actual = moments( edges_current_graph( inneriter,1 ), edges_current_graph( inneriter,2 ) );
                    m_projection = moments_edges_removed( inneriter );
                    curr_edge_decrease_value = 0.5 * ( ( ( 1 + m_actual ) * log( ( 1 + m_actual ) / ( 1 + m_projection ) ) ) + ( ( 1 - m_actual ) * log( ( 1 - m_actual ) / ( 1 - m_projection ) ) ) );
                    curraicimprovement = edge_importance_scaledownfactor - curr_edge_decrease_value;
                    if(curraicimprovement > maxaicimprovement)
                        maxaicimprovement = curraicimprovement;
                        maxaicimprovementedge = [ edges_current_graph( inneriter,1 ), edges_current_graph( inneriter,2 ) ];
                    end
                end
            elseif( iter > 2 && pruneonthego )
                for inneriter = 1 : size(edges_current_graph,1)
                    if( ( previous_added_edge(1) == edges_current_graph(inneriter,1) && ...
                            previous_added_edge(2) == edges_current_graph(inneriter,2) ) || ...
                        ( previous_added_edge(2) == edges_current_graph(inneriter,1) && ...
                            previous_added_edge(1) == edges_current_graph(inneriter,2) ) )
                        thresholdthetacurrround = theta_current_graph(inneriter);
                        break;
                    end
                end

                for inneriter = 1 : length(theta_current_graph)
                    if(abs(theta_current_graph(inneriter)) < edge_importance_scaledownfactor * abs(thresholdthetacurrround))
                        undirected_graph(edges_current_graph(inneriter,1),edges_current_graph(inneriter,2)) = 0;
                        undirected_graph(edges_current_graph(inneriter,2),edges_current_graph(inneriter,1)) = 0;
                    end
                end
                edge_importance_scaledownfactor = edge_importance_scaledownfactor * ...
                    edge_importance_scaledownfactor_constant;

                [moments_current_graph,edges_current_graph,theta_current_graph] = get_projection( undirected_graph, moments );

                remaining_edges = tril( complete_graph ) - tril(undirected_graph);
            end
            
            [ rows_remaining_edges, cols_remaining_edges ] = find( remaining_edges );

            edge_improvement = zeros(length(rows_remaining_edges), 3);
            for index = 1 : length( rows_remaining_edges )
                m_actual = moments( rows_remaining_edges( index ), cols_remaining_edges( index ) );
                m_projection = moments_current_graph( rows_remaining_edges( index ), cols_remaining_edges( index ) );
                curr_edge_improvement_value = 0.5 * ( ( ( 1 + m_actual ) * log( ( 1 + m_actual ) / ( 1 + m_projection ) ) ) + ( ( 1 - m_actual ) * log( ( 1 - m_actual ) / ( 1 - m_projection ) ) ) );
                if( rows_remaining_edges(index) == 1 || cols_remaining_edges(index) == 1 )
                    curr_edge_improvement_value = meanweight * curr_edge_improvement_value;
                end
                curr_edge_improvement = [ rows_remaining_edges( index ), cols_remaining_edges( index ), curr_edge_improvement_value ];
                %edge_improvement = [ edge_improvement; curr_edge_improvement ];
                edge_improvement(index, :) = curr_edge_improvement;
            end

            %Sort in descending order (and hence the -)
            edge_improvement = sortrows( edge_improvement, -3 );
            %valid_edge_improvement = ~isnan(edge_improvement(:,3));
%             algo_progress_profile{iter} = [ edges_already_added; edge_improvement(valid_edge_improvement,:) ];
            theta_profile{iter} = theta_current_graph;
            edges_profile{iter} = edges_current_graph;

            edge_index = 1;
            while( edge_index <= size(edge_improvement,1) )
                if ~isnan(edge_improvement( edge_index,3 ))
                    break;
                end
                remaining_edges( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 0;
                edge_index = edge_index+1;
            end
            %If the highest improvement is very less, do not add any more edges
            if( edge_index > size(edge_improvement,1) || edge_improvement(edge_index,3) < improvement_threshold )
                break;
            end

            while( edge_index <= size( edge_improvement, 1 ) )
                remaining_edges( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 0;
                undirected_graph( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 1;
                undirected_graph( edge_improvement( edge_index, 2 ), edge_improvement( edge_index, 1 ) ) = 1;
                % We assume that if adding a particular edge makes the graph
                % non-planar, then the corresponding moment will be 'NaN'
                % If not, we need to do a planarity test here again
                if( isnan( edge_improvement( edge_index, 3 ) ) )
                    undirected_graph( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 0;
                    undirected_graph( edge_improvement( edge_index, 2 ), edge_improvement( edge_index, 1 ) ) = 0;
                elseif( addpruneonthego )
                    if( ( edge_improvement(edge_index,3) - edge_importance_scaledownfactor > 0 ...
                            || maxaicimprovement > 0 ) )
                        if(edge_improvement(edge_index,3) - edge_importance_scaledownfactor >= maxaicimprovement )
                            graph( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 1;
                            %edges_already_added = [ edges_already_added ; ...
                            %            edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ), -1 ];
                            %edge_addition_order = [ edge_addition_order; edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ];
                            edge_addition_order(iter, :) = [edge_improvement(edge_index, 1), edge_improvement(edge_index, 2)];
                            nllimprovement_lowerbound_info(iter) = edge_improvement(edge_index,3);
                            edges_info(iter,:) = [ edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ), 1 ];
                            break;
                        else
                            undirected_graph( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 0;
                            undirected_graph( edge_improvement( edge_index, 2 ), edge_improvement( edge_index, 1 ) ) = 0;
                            undirected_graph(maxaicimprovementedge(1),maxaicimprovementedge(2)) = 0;
                            undirected_graph(maxaicimprovementedge(2),maxaicimprovementedge(1)) = 0;
                            remaining_edges = tril( complete_graph ) - tril(undirected_graph);
                            nllimprovement_lowerbound_info(iter) = maxaicimprovement - edge_importance_scaledownfactor;
                            edges_info(iter,:) = [ maxaicimprovementedge(1), maxaicimprovementedge(2), -1 ];
                            break;
                        end
                    else
                        remaining_edges( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 1;
                        undirected_graph( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 0;
                        undirected_graph( edge_improvement( edge_index, 2 ), edge_improvement( edge_index, 1 ) ) = 0;

                        edge_importance_scaledownfactor = edge_importance_scaledownfactor * ...
                            edge_importance_scaledownfactor_constant;
                        break;
                    end
                else
                    graph( edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ) = 1;
                    %edges_already_added = [ edges_already_added ; ...
                    %            edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ), -1 ];
                    %edge_addition_order = [ edge_addition_order; edge_improvement( edge_index, 1 ), edge_improvement( edge_index, 2 ) ];
                    edge_addition_order(iter, :) = [edge_improvement(edge_index, 1), edge_improvement(edge_index, 2)];
                    previous_added_edge = [ edge_improvement(edge_index,1), edge_improvement(edge_index,2) ];
                    if( edge_improvement(edge_index,3) < lowest_improvement_threshold )
                        lowest_improvement_threshold = edge_improvement(edge_index,3);
                    end
                    break;
                end
                edge_index = edge_index+1;
            end
        end
        
        algo_progress_profile = {edges_info; nll_info(2:end); nllimprovement_lowerbound_info};
        %algo_progress_profile{1} = edges_info;
        %algo_progress_profile{2} = nll_info;
        %algo_progress_profile{3} = nllimprovement_lowerbound_info;
%         continueiterations = false;
%         for count = 1 : size(edges_current_graph,1)
%             current_theta = theta_current_graph(count);
%             current_edge_final_benefit = current_theta * tanh(current_theta) - log(cosh(current_theta));
%             if( current_edge_final_benefit < lowest_improvement_threshold )
%                 continueiterations = true;
%                 undirected_graph(edges_current_graph(count,1),edges_current_graph(count,2)) = 0;
%                 undirected_graph(edges_current_graph(count,2),edges_current_graph(count,1)) = 0;
%             end
%         end
%         graph = tril( undirected_graph );
%         complete_graph = ones( n ) - eye( n );
%         remaining_edges = tril( complete_graph ) - graph;
%     end