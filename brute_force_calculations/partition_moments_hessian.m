function [Z,M,H] = partition_moments_hessian(n,edges,theta)

    Z = 0;
    M_edges = zeros( 1, size(edges,1) );
    M = zeros( n );
    H = zeros( size(edges,1) );
    for index = 0 : 2^n - 1
        curr_Z = 1;
        num = dec2bin(index,n);
        for count = 1 : size(edges,1)
            if( num(edges(count,1)) == num(edges(count,2)) )
                curr_Z = curr_Z * exp(theta(count));
            else
                curr_Z = curr_Z * exp(-theta(count));
            end
        end
        Z = Z + curr_Z;
        
        for i = 1 : size(edges,1)
            M_edges(i) = M_edges(i) + (2*bin2dec(num(edges(i,1)))-1) * (2*bin2dec(num(edges(i,2)))-1) * curr_Z;
        end
        
        for i = 1 : n
            for j = 1 : n
                M(i,j) = M(i,j) + (2*bin2dec(num(i))-1) * (2*bin2dec(num(j))-1) * curr_Z;
            end
        end
        
        for i = 1 : size(edges,1)
            for j = i+1 : size(edges,1)
                H(i,j) = H(i,j) + (2*bin2dec(num(edges(i,1)))-1) * (2*bin2dec(num(edges(i,2)))-1) * ...
                                  (2*bin2dec(num(edges(j,1)))-1) * (2*bin2dec(num(edges(j,2)))-1) * ...
                                  curr_Z;
            end
        end
    end
    
    M = M / Z - eye(n);
    M_edges = M_edges / Z;
    H = eye(size(edges,1)) + ( ( H + H')/Z ) - M_edges'*M_edges;