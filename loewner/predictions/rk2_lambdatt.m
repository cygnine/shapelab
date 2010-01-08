function[lambdatt] = rk2_lambdatt(a, g, ga, ds, lambdat)
% rk2_lambdatt -- The relation connecting lambdat and lambdatt
%
% ltt = rk2_lambdatt(a, g, lambdat)

lambdatt = (4*imag(a)/ds - imag(g)*lambdat)./(3*imag(ga));
