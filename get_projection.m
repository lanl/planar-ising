% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function [moments_current_graph,edges,theta] = get_projection( graph, moments )

%Back-tracking line search parameters
alpha = 0.4;
eta = 0.3;
threshold = 7e-5;%7e-4;
threshold_stop_backtracking = 1e-1;
deltatheta_max_threshold = 5;%2;
%Truncate the moments to those on the edges of the graph
n = length( graph );
moments_truncated = [];
edges = [];
for j = 1 : n
    for i = j+1 : n
        if( graph(i,j) == 1)
            edges = [ edges ; [ i, j ] ];
            moments_truncated = [ moments_truncated, moments(i,j) ];
        end
    end
end
theta = zeros( 1, size(edges, 1) );

node_positions = chrobak_payne_straight_line_drawing( graph );

% Implement gradient descent method to get to the minimum of the convex fn :
% log( Z( theta ) ) - mu.theta
iter = 0;
while true
    iter = iter+1;
	[ partition_function, moments_current_graph, hessian ] = kacward( node_positions, edges, theta, true );
    
    if( abs(cond(hessian)) > 1e10 )
        error( 'hessian is nearly singular' );
    end

	descent_step = 1;
	moments_diff = moments_current_graph - moments_truncated;
    %delta_theta = - ( descent_step * moments_diff * inv(hessian)' );
    delta_theta = - descent_step * (hessian\moments_diff')';
    if( abs(sum( moments_diff .* delta_theta )) < threshold )
		break;
    elseif( abs(sum( moments_diff .* delta_theta )) > threshold_stop_backtracking )
        % Decrease the stepping size till the difference between the function values before and after the jump is atleast
        % alpha times what the slope suggests it will be
        while true
            if( sum(abs(delta_theta)) < 1e-30 )
                error('this should not happen');
                break;
            end
            if( max(abs(delta_theta)) > deltatheta_max_threshold )
                delta_theta = delta_theta / max(abs(delta_theta));
            end
            theta_next_guess = theta + delta_theta;
            [ partition_function_next_step, ~ ] = kacward( node_positions, edges, theta_next_guess, false );
            if( negative_likelihood( partition_function_next_step, moments_truncated, theta_next_guess ) < negative_likelihood( partition_function, moments_truncated, theta ) + alpha * descent_step * sum( moments_diff .* delta_theta ) )
                break;
            end
    		%descent_step = descent_step * eta;
            %delta_theta = - ( descent_step * moments_diff * inv(hessian)' );
            %delta_theta = - descent_step * (hessian\moments_diff')';
            delta_theta = eta * delta_theta;
        end
    end

	%theta = theta - ( descent_step * moments_diff * inv(hessian)' );
    theta = theta + delta_theta;
    
end

% zero_theta_edges = find( abs(theta) < exp(-10) );
% for index = 1 : length( zero_theta_edges )
%     a( edges( zero_theta_edges( index ), 1 ), edges( zero_theta_edges( index ), 2 ) ) = 0;
%     a( edges( zero_theta_edges( index ), 2 ), edges( zero_theta_edges( index ), 1 ) ) = 0;
% end
% 
% non_zero_theta_edges = find( abs(theta) > exp(-10) );
% theta = theta( non_zero_theta_edges );
% edges = edges( non_zero_theta_edges, : );

[~, moments_current_graph] = kacward_allmoments( graph, edges, theta );




