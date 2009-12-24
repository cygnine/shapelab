function[s,lambda,dlambda] = compute_driving_function(z, varargin)
% compute_driving_function -- point-evaluations of a Loewner ODE driving function
%
% [s,lambda,dlambda] = compute_driving_function(z, {s=[]},visualize=false)
%
%     Given N samples z \in \mathbbm{H}, this function returns a length-N vector
%     with point-evaluations of the driving function lambda(s). The
%     parameterization s is determined on-the-fly and depends on the locations
%     of the points z. We assume z(1) is real -- that it is on the boundary of
%     \mathbbm{H}.
%
%     If we view the samples z as coming from the trace of a curve \gamma, then
%     the driving function can be used to evolve the Loewner ODE to generate a
%     map from \mathbbm{H} without \gamma onto \mathbbm{H}. I.e., it solves a
%     conformal mapping problem.

persistent strict_inputs invert_a
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.loewner import invert_a
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
g = zeros([N 1]);
lambda = zeros([N 1]);
dlambda = zeros([N-1 1]);
branch = true([N 1]);
s = zeros([N 1]);

% Initial data
g = z;
lambda(1) = z(1);
a = 1/2*g.^2 - lambda(1)*g;

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
  % Recover g:
  g = invert_a(a, lambda(q+1),branch);

  if opt.visualize
    set(myplot, 'xdata', real(g), 'ydata', imag(g));
    drawnow;
  end
end
