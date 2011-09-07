function[thetas] = epdiff_particles(theta, a, b, ts, varargin)
% epdiff_particles -- Evolves particle positions using epdiff
%
% thetas = epdiff_particles(theta, a, b, ts, {dt=0.01, wrap=false, c=zeros([3 1])})
%
%     Positions initial teichons at locations b \in S^1 and strengths a. Then
%     epdiff is used to evolve the teichon positions. However, the initial
%     locations (particles) theta are affected the by the velocity field induced
%     by the teichons. This function evolves the positions of the particles
%     under the (self-updating) velocity field produced by the teichons.
%     Therefore, this function implicitly solves epdiff in order to evolve the
%     teichons.
%
%     The output thetas is a matrix, where column k corresponds to the particle
%     positions at time ts(k).
%
%     Optional inputs dt (for time-stepping of epdiff) and wrap are given to
%     indicate whether the values of the particles on S^1 are `wrapped' to [0,
%     2*pi) or not. The final optional input is a 3-vector that corresponds to a
%     coordinates of the element of the velocity field in the kernel of the WP
%     operator L. The functional form of the kernel basis is dictated by
%     shapelab.wp.kernel_basis.

persistent strict_inputs teichon_wp epdiff_teichon_rhs
persistent rk interval_wrap wrap kernel_basis greens_function
if isempty(rk)
  from labtools import strict_inputs 
  from labtools import interval_wrap 
  from shapelab.wp import teichon_wp epdiff_teichon_rhs kernel_basis greens_function

  from odesolve.coeffs import lserk4
  rk = lserk4();
  wrap = @(t) interval_wrap(t, [0, 2*pi]);
end

N = length(a);
if N < 4
  fprintf('With less than four samples, the evolution is trivial...aborting\n');
end

opt = strict_inputs({'dt', 'wrap', 'c'}, {0.01, false, zeros([3 1])}, [], varargin{:});

theta = theta(:); 
T = max(ts);

dt = opt.dt;
t = 0:dt:T;
if abs(T-t(end))>1e-13
  t(end+1) = T;
end

t = sort([t(:); ts(:)]);
flags = find(diff(t)<1e-12);
t(flags+1) = [];
% Vector of teichon information:
y = zeros([2*N 1]);

% Initial data
y(N+1:end) = b;
y(1:N) = a;
if opt.wrap
  y(N+1:end) = wrap(y(N+1:end));
end

time = 0;
step_number = 1;
flags = abs(time-ts)<1e-12;
if any(flags)
  thetas(:,flags) = theta;
end

% Implementation:
while time<T;
  k = zeros([2*N 1]);
  ktheta = zeros(size(theta));
  stage_time = t;
  dt = t(step_number+1) - t(step_number);

  % RK stages
  for p = 1:length(rk.a);
    stage_time = time + dt*rk.c(p);

    f = epdiff_teichon_rhs(y(1:N),y(N+1:end));
    [t1,t2] = meshgrid(y(N+1:end), theta);
    ftheta = greens_function(t1-t2)*y(1:N) + kernel_basis(theta)*opt.c;

    k = rk.a(p)*k + dt*f;
    ktheta = rk.a(p)*ktheta + dt*ftheta;

    y = y + rk.b(p)*k;
    theta = theta + rk.b(p)*ktheta;

    if opt.wrap
      %y(N+1:end) = wrap(y(N+1:end));
      theta = wrap(theta);
    end

  end

  % Update time, iteration count
  time = t(step_number+1);
  step_number = step_number+1;

  % Give output
  flags = abs(time-ts)<1e-12;
  if any(flags)
    thetas(:,flags) = theta;
  end
end
