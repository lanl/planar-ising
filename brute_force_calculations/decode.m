% Jason K. Johnson, MIT, Oct 2006.

function A = decode(a)

% SYNOPSIS A = decode(a)
% Decode 32-bit unsigned integer as a subset of {0,...,31}.

A = find(bitget(a,1:32))-1;
