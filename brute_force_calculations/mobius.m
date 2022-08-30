% Jason K. Johnson, MIT, Oct 2006.

function g = mobius(omega,f)

N = length(f);
n = round(log2(N));
if (2^n ~= N)
  error('length of f must be power of two');
end

g = f;
for k=0:n-1
  g = reshape(g,2^k,N/(2^k));
  g(:,2:2:end) = g(:,2:2:end) + omega * g(:,1:2:end);
end
g = reshape(g,size(f));
