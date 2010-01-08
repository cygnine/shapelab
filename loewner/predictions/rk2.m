function[data] = rk2(data)
% rk2 -- Loewner evolution predictions using an RK2 approximation
%
% data = rk2(data)
%
%     Predicts the value of the derivative of the driving function using an RK2
%     evolution approximation. This prediction comes in the form of the
%     derivative of lambda at time-points 0*ds and 2/3*ds.

persistent bisection
if isempty(ga)
  %from shapelab.loewner.predictions import rk2_g_aprime as ga
  %from shapelab.loewner.predictions import lambdatt
  from labtools.rootfind import bisection 
end

% Subfucntion
function[y] = f(x, a, g, ds, lambda)

persistent lambdatt ga
if isempty(lambdatt)
  from shapelab.loewner.predictions import rk2_lambdatt as lambdatt
  from shapelab.loewner.predictions import rk2_g_aprime as ga
end

gap = ga(a, g, ds, lambda, x);
lambdattx = lambdatt(a, g, gap, ds, x);
lambdanp1 = lambda + ds/4*(x + 3*lambdattx);

y = real(a) + ds/4*(2 - x*.real(g) + 6 - 3*lambdattx.*real(gap)) + 1/2*lambdanp1.^2;
