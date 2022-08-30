% Jason K. Johnson, MIT, Oct 2006.

function theta = invert_moments(eta)

p_eta = mobiusT(-1,eta);
theta = mobius(-1,log(p_eta));


