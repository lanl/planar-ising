% Jason K. Johnson, MIT, Oct 2006.

function a = encode(A)
% SYNOPSIS a = encode(A)
% Encodes a subset of 0,...,31 as a 32-bit unsigned integer.

a = uint32(sum(uint32(2.^uint32(A))));
