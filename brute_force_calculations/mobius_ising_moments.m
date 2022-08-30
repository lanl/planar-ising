function [ moments_nodes_ising, moments_nodepairs_ising ] = mobius_ising_moments(number_of_nodes,edges,theta_nodes,theta_edges)

    number_of_edges = size(edges,1);
    [beta_nodes,beta_edges,beta_const] = ising2boltz(edges,theta_nodes,theta_edges);

    theta_vector = zeros(1,2^(number_of_nodes));

    theta_vector(1+encode([])) = beta_const;

    for i = 1 : number_of_nodes
        theta_vector(1+encode([i-1])) = beta_nodes(i);
    end

    for i = 1 : number_of_edges
        theta_vector(1+encode(edges(i,:)-1)) = beta_edges(i);
    end

    moments_boltzman = moments( theta_vector );

    moments_nodes_ising = zeros(1,number_of_nodes);
    for i = 1 : number_of_nodes
        moments_nodes_ising(1,i) = moments_boltzman(1+encode([i-1]));
    end
    
    moments_nodepairs_ising = zeros(number_of_nodes);
    for i = 1 : number_of_nodes
        for j = i+1 : number_of_nodes
            moments_nodepairs_ising(i,j) = moments_boltzman(1+encode([i-1,j-1]));
        end
    end

    %moments_nodepairs_ising = moments_nodepairs_ising + moments_nodepairs_ising';
    %moments_nodepairs_ising = 4 * moments_nodepairs_ising - ( ones(number_of_nodes) - eye(number_of_nodes) );
    
    for i = 1 : number_of_nodes
        for j = i+1 : number_of_nodes
            moments_nodepairs_ising(i,j) = 1+4*moments_nodepairs_ising(i,j)...
                                            -2*moments_nodes_ising(1,i)-2*moments_nodes_ising(1,j);
        end
    end
    moments_nodepairs_ising = moments_nodepairs_ising + moments_nodepairs_ising';
    
    moments_nodes_ising = 2*moments_nodes_ising - 1;
