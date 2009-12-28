function[s,lambda,dlambda,g,gn] = compute_driving_unzip(z, varargin)
% compute_driving_unzip -- point-evaluations of a Loewner ODE driving function
%
% [s,lambda,dlambda,gplus,gminus] = compute_driving_unzip(z, {visualize=false})
%
%     This function performs the same operations as compute_driving_function,
%     except it also gives the outputs gplus and gminus, which are evaluations
%     of the terminal conformal map of the locations z. This is used to obtain
%     the welding map for shapes.

persistent strict_inputs invert_a
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.loewner import invert_a_unzip as invert_a
end

z = z(:);
N = length(z);

opt = strict_inputs({'visualize'}, {false}, [], varargin{:});
assert(imag(z(1))<eps, 'Error: the first sample point must be real');

% The differential equation is:
%  a' = 2 - g*lambda'
%  lambda' = lambda'
%
% The variable g now represents the evaluation of the conformal map

% Storage allocation:
a = zeros([N 1]);
an = zeros([N 1]);
g = zeros([N 1]);
gn = zeros([N 1]);
lambda = zeros([N 1]);
dlambda = zeros([N-1 1]);
s = zeros([N 1]);

% Initial data
g = z;
gn = g;
lambda(1) = z(1);
a = 1/2*g.^2 - lambda(1)*g;
an = a;

if opt.visualize
  myplot = plot(g, 'r.');
end

% For now, use Euler-stepping
for q = 1:(N-1)
  % Solution of a linear system enforcing
  %  (1): real-valued g(q+1)
  %  (2): coincidence of g(q+1) and lambda(q+1)
  dl = imag(a(q+1))/imag(g(q+1));
  ds = real((dl^2 + 2*lambda(q)*dl - 2*g(q+1)*dl + lambda(q)^2 + g(q+1)^2 - ...
             2*lambda(q)*g(q+1))/(-4));

  s(q+1) = s(q) + ds;
  
  % An approximation to derivative of lambda:
  dlambda(q) = dl/ds;

  % Evolve lambda:
  lambda(q+1) = lambda(q) + dl;

  % Evolve the rest of the points:
  a = a + ds*(2-g.*dlambda(q));

  a(q+1) = -1/2*(lambda(q+1).^2);  % Explicitly enforce -- machine roundoff crap

  an(q+1:end) = a(q+1:end);
  an(1:q) = an(1:q) + ds*(2-gn(1:q).*dlambda(q));

  % Recover g:
  [g,gn] = invert_a(a, an, lambda(q+1));

  if opt.visualize
    set(myplot, 'xdata', real(g), 'ydata', imag(g));
    drawnow;
  end
end
