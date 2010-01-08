function[g] = rk2_g_aprime(a, g, ds, lambda, dlambda)
% rk2_g_aprime -- An approximation to g(a') via evolution
%
% g = rk2_g_aprime(a, g, dlambda)

persistent invert_a
if isempty(invert_a)
  from shapelab.loewner import invert_a
end

ap = a + 2/3*ds*(2-g*dlambda);

lambda_temp = lambda + 2/3*ds*dlambda;
g = invert_a(ap, lambda_temp, true(size(dlambda)));
