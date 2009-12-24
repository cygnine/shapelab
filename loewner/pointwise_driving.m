function[g] = pointwise_driving(s, lambda, dlambda, z, varargin)
% pointwise_driving -- Evolution of the Loewner DE using point-evaluations
%
% g = pointwise_driving(s, lambda, dlambda, z, {visualize=false})
%
%   Given a collection of sample points z, this function evolves the Loewner
%   equation using the parameterization lambda(s) that is given as input. If s
%   and lambda are of length N, then dlambda is an N-1 length vector with
%   evaluations of the derivative of lambda at locations s(1:N-1). The output is
%   g = g(z), evaluations of the conformal map taking a slit approximation
%
%   As of now, this function just uses the forward Euler approximation.

persistent strict_inputs invert_a
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.loewner import invert_a
end

opt = strict_inputs({'visualize'}, {false}, [], varargin{:});

zsize = size(z);
z = z(:);
N = length(z);

% Allocation:
a = zeros([N 1]);
g = zeros([N 1]);
branch = true([N 1]);

% Initial data
g = z;
a = 1/2*g.^2 - lambda(1)*g;

ds = diff(s);

if opt.visualize
  myplot = plot(g, 'r.');
end

for q = 1:length(ds)
  a = a + ds(q)*(2-g.*dlambda(q));
  g = invert_a(a, lambda(q+1), branch);

  if opt.visualize
    set(myplot, 'xdata', real(g), 'ydata', imag(g));
    drawnow
  end
end

g = reshape(g, zsize);
