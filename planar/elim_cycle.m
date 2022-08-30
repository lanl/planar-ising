

n = 4;

K = zeros(n,2*n);

for i = 1:n
  K(i,2*(i-1)+1) = 1;
  K(i,2*(i-1)+2) = 1;
end

K

C = zeros(2*n,2*n); 

for i = 1:(2*n-1)
  C(i,i+1) = 2*mod(i+1,2)-1;
end

C(1,2*n) = 2*mod(n+1,2)-1;

C = C - C';

C

Z_C = det(C)

A_elim = K*inv(C)*K'
