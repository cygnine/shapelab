function[cint,sint] = trig_integral2(k,s)
% trig_integral2 -- evaluates a Fresnel-like trigonometric integral
%
% [cint,sint] = trig_integral2(k,s)
%
% For a 2-vector s, with k of length N, this function evaluates 
%
% cint(1) = \int_s(1)^s(2) cos(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
% cint(2) = \int_s(1)^s(2) s*cos(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
% cint(3) = \int_s(1)^s(2) s^2*cos(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
%    ...
% cint(N) = \int_s(1)^s(2) s^(N-1)*cos(k(1) + k(2)*s + k(3)*s^2 + ... ) ds

% sint(1) = \int_s(1)^s(2) sin(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
% sint(2) = \int_s(1)^s(2) s*sin(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
% sint(3) = \int_s(1)^s(2) s^2*sin(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
%    ...
% sint(N) = \int_s(1)^s(2) s^(N-1)*sin(k(1) + k(2)*s + k(3)*s^2 + ... ) ds
%
% This function is vectorized in the rows of s, so s can be an M x 2 vector, and
% then the outputs cint and sint are size (M x N) matrices.

persistent x w gq peval repnodes N spdiag
if isempty(gq)
  from speclab.orthopoly1d.jacobi.quad import gauss_quadrature as gq
  from speclab.monomials import evaluate as peval
  from piecewise_interpolation.grid_tools import replicate_local_nodes as repnodes
  from labtools import spdiag

  % Default order of quadrature: 
  N = 50;
  [x,w] = gq(N, 'alpha', 0, 'beta', 0);
end

K = size(s,1);
order = length(k);
scale = diff(s,1,2).'/2;
shift = mean(s,2).';

xs = repnodes(x, s);

%xs = repmat(x, [1 K]);
%xs = xs*spdiags(scale(:),0,K,K) + repmat(shift, [N 1]);

xs = x*ones([1 K]);
xs = xs*spdiag(scale(:)) + ones([N 1])*shift;

k = flipud(k(:));

p = polyval(k, xs);

cint = zeros([K order]);
sint = zeros([K order]);

temp = spalloc(N,N,N);
inds = (1:N) + (0:(N-1))*N;
temp(inds) = w;

cp = temp*cos(p);
sp = temp*sin(p);

for q = 1:order
  cint(:,q) = (sum(cp,1).*scale).';
  sint(:,q) = (sum(sp,1).*scale).';

  cp = cp.*xs;
  sp = sp.*xs;
end
