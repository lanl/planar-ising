% Jason K. Johnson, MIT, Oct 2006.

function eta = moments(theta)

p_theta = exp(mobius(1,theta));
eta = mobiusT(1,p_theta);
eta = eta / sum(p_theta);
