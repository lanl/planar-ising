% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

function [a, g, node_positions] = outerplanargraph()
n = 12;
a = zeros(n);
for i = 1 : n-1
a(i,i+1)=1;
end
a(1,n)=1;
a(1,3)=1;a(1,4)=1;a(1,11)=1;a(4,11)=1;a(5,11)=1;a(6,11)=1;a(6,8)=1;a(6,9)=1;a(9,11)=1;

a = a+a';

angles = [ -pi : 2*pi/n : (n-2)*pi/n ]';
node_positions = n * [ cos(angles), sin(angles) ];

% need weights
tmp = -0.95 + 1.9 * rand(n);
tmp = tmp + sign(tmp)*0.05;
tmp = triu(tmp,1);
tmp = tmp + tmp';
g = a .* tmp;
clear tmp
