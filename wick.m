% ©2022. Triad National Security, LLC. All rights reserved.
% This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.



function K = wick(n,edges,mu,H)

K = zeros(n,n);
Gamma = zeros(n,n);

s = edges(:,1);
t = edges(:,2);

arcs = [[s t]; [t s]];

mu = [mu; mu];
H = [H H; H H];

for i = 1:n
  K(i,i) = 1;
  Gamma(i,i) = 1;
end

m = length(arcs);

for ab = 1:m
  a = arcs(ab,1);
  b = arcs(ab,2);
  K(a,b) = mu(ab);
  Gamma(a,b) = 1;
end

while (any(any(Gamma==0)))

  for ab = 1:m
    for cd = 1:m

      a = arcs(ab,1);
      b = arcs(ab,2);
      c = arcs(cd,1);
      d = arcs(cd,2);
 
      if (~Gamma(a,d) & Gamma(a,c) & Gamma(b,d) & Gamma(b,c))
        K(a,d) = (H(ab,cd)-K(a,c)*K(b,d))/K(b,c);
        Gamma(a,d) = 1;
      end
 
    end
  end

end