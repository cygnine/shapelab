function[zc] = product_driving(xi, T, z, varargin)
% product_driving -- Backwards solution to the Loewner equation
%
% zc = product_driving(xi, T, {dt=0.01, visualize=false})
%
%     Solves the Loewner evolution equation
%
%       dg/dt = 2/(g - xi),
%
%     backwards for every z in the input vector z, where g is a conformal map from H \
%     (gamma) to H. xi is a function handle pointing to a function of one
%     variable, t. It is the "driving function".  The output zc is a collection of
%     evaluations of g^{-1}([0,T]); i.e. it is the preimage of the range of the
%     driving function.

persistent input_schema lserk4 eno_reconstruction
if isempty(input_schema)
  from labtools import input_schema
  from odesolve.coeffs import lserk4
  from eno import eno_reconstruction
end

rk = lserk4();

inputs = {'dt', 'visualize'};
defaults = {0.01, false};

opt = input_schema(inputs, defaults, [], varargin{:});

% Find the derivative of xi
x = linspace(0, T, 100);
y = xi(x);
dxi = diff(eno_reconstruction(x, y, 'k', 5));

% Reflect xi so rk can work forward in time
xit = @(t) xi(T - t);
t = 0;
niter = 0;

N = ceil(T/opt.dt);
zsize = size(z);
z = z(:);
zc = z;
% Initial value of a:
%a = zc.^2 - 2*zc*xit(0);
a = 1/2*zc.^2 - zc*xit(0);

if opt.visualize
  plot_id = make_plots();
  update_plots()
end

while t < T
  ku = zeros(size(a));

  if t+opt.dt>T;
    opt.dt = T - t; % to evolve exactly to T.
  end

  ts = t + rk.c*opt.dt;
  dxis = dxi(ts);

  % RK stages
  for p = 1:length(rk.a)
    %a_rhs = -2*(1 + zc.*dxis(p));
    a_rhs = -2 + zc.*dxis(p);

    ku = rk.a(p)*ku + opt.dt*a_rhs;

    a = a + rk.b(p)*ku;
    if niter==0;
      zc = invert_a(a, zc, xit(ts(p)),'clear');
    else
      zc = invert_a(a, zc, xit(ts(p)));
    end

  end

  if opt.visualize
    update_plots()
  end

  t = t+opt.dt;
  niter = niter+1;

end

zc = reshape(zc, zsize);

%%%%%%%%%%%%%% Nested functions %%%%%%%%%%%%%%%%%

function plot_id = make_plots()
  plot_id = plot(z, 'r.');
end

function [] = update_plots()
  set(plot_id, 'xdata', real(zc), 'ydata', imag(zc));
  axis([-5, 5, 0, 5]);
  drawnow; pause(0.01)
end

end

%%%%%%%%%%%%%% Subfunctions %%%%%%%%%%%%%%%%%

function z = invert_a(a, g, xi, varargin);
% invert_a -- inverts the a(g) relation
%
% z = invert_a(a, g, xi, {'clear'})
%
%     The relation 2*a = z^2 - 2*z*xi is inverted to solve for z, with appropriate
%     choices of branches for z. The input g is the previous value of z:
%     necessary to choose the correct branch.
%
%     The optional clear argument tells the routine to clear previously saved
%     data.

persistent csqrt cflags pflags remainder discriminant

discriminant = 2*a + xi.^2;
% discriminant = a + xi.^2;

if isempty(csqrt) | not(isempty(varargin))
  from shapelab.common import positive_angle_exponential as csqrt
  cflags = false(size(a));
  pflags = false(size(a));
  remainder = false(size(a));

  pflags = abs(imag(discriminant))<1e-14;
  cflags = not(pflags);
end

z = zeros(size(a));

temp = discriminant(pflags)>0;
remainder(pflags) = not(temp);
pflags(pflags) = discriminant(pflags)>0;

z(cflags) = xi + csqrt(discriminant(cflags), 1/2);

signs = sign(g(pflags) - xi);
z(pflags) = xi + signs.*sqrt(discriminant(pflags));

z(remainder) = xi + i*sqrt(abs(discriminant(remainder)));

cflags(remainder) = true;
remainder(:) = false;

end

