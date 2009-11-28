function[cint,sint] = trig_integral(k,s)
% trig_integral -- evaluates a Fresnel-like trigonometric integral
%
% [cint,sint] = trig_integral(k,s)
%
% For a 2-vector s, this function evaluates 
%
% cint = \int_s(1)^s(2) cos(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
% sint = \int_s(1)^s(2) sin(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
%
% This function is vectorized in the rows of s, so s can be an N x 2 vector, and
% then the outputs c and s are length-N vectors.

persistent x w gq peval repnodes N
if isempty(gq)
  from speclab.orthopoly1d.jacobi.quad import gauss_quadrature as gq
  from speclab.monomials import evaluate as peval
  from piecewise_interpolation.grid_tools import replicate_local_nodes as repnodes

  % Default order of quadrature: 
  N = 100;
  [x,w] = gq(N, 'alpha', 0, 'beta', 0);
end

K = size(s,1);
scale = diff(s,1,2).'/2;
shift = mean(s,2).';

xs = repnodes(x, s);
xs = repmat(x, [1 K]);
xs = xs*spdiags(scale(:),0,K,K) + repmat(shift, [N 1]);

k = flipud(k(:));

p = polyval(k, xs);

cint = ((w.'*cos(p)).*scale).';
sint = ((w.'*sin(p)).*scale).';
