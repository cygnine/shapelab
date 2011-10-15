function[a_out, b_out, y_out, wp_values] = epdiff_teichon(a0, b0, y0, ts, varargin)
% epdiff_teichon -- Evolves epdiff using a teichon ansatz
%
% [a_out, b_out, y_out, wp_values] = epdiff_teichon(a0, b0, y0, ts, {dt=0.01,
%                                                   wrap=false})
%
%     Solves epdiff on elements in Diff(S^1) using a teichon ansatz. The inputs
%     a0 and b0 are the initial strengths and locations of the teichons; y0 are
%     the initial locations of particles that one wishes to track (may be
%     empty).  The teichons have velocity representation:
%
%      v(theta) = sum_j a_j G(theta-b_j) + sum_q k_q phi_q(theta),
%
%     where G is the WP-Green's function, and phi_q for q = [1,2,3] is a basis
%     for the kernel of the WP norm (see shapelab.wp.wp_minimal_metric and
%     shapelab.wp.kernel_basis). 
%
%     Evolution proceeds using an ODE for the positions and locations of the
%     teichons; a 5-stage 4th order SSP LSERK method is used. The outputs are
%     returned at the t values in the input array ts. The initial time is
%     assumed to be t=0.
%
%     If the optional input wrap is set to true, then the locations b of the
%     teichons are wrapped to the [0, 2*pi) interval.

persistent input_parser parser
persistent teichon_wp epdiff_teichon_rhs
persistent rk interval_wrap wrap
if isempty(parser)
  from labtools import input_parser 
  from labtools import interval_wrap 
  from shapelab.wp import teichon_wp epdiff_teichon_rhs

  [opt,parser] = input_parser({'dt', 'wrap'}, {0.01, false}, [], varargin{:});

  from odesolve.coeffs import lserk4
  rk = lserk4();
  wrap = @(t) interval_wrap(t, [0, 2*pi]);
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

N = length(a0);
M = length(y0);
if N < 4
  fprintf('With less than four samples, the evolution is trivial...aborting\n');
  a_out = [];
  b_out = [];
  y_out = [];
  k_out = [];
  wp_values = [];
end

T = max(ts);

% Randomly: dt = 0.01
dt = opt.dt;
t = 0:dt:T;
if abs(T-t(end))>1e-13
  t(end+1) = T;
end

t = sort([t(:); ts(:)]);
flags = find(diff(t)<1e-12);
t(flags+1) = [];
y = zeros([2*N+M 1]);

% Initial data
y(1:N) = a0;
y(N+1:2*N) = b0;
y(2*N+1:end) = y0;
if opt.wrap
  y(N+1:end) = wrap(y(N+1:end));
end

time = 0;
step_number = 1;
wp_values = zeros(size(ts));
flags = abs(time-ts)<1e-12;
if any(flags)
  a_out(:,flags) = y(1:N);
  b_out(:,flags) = y(N+1:2*N);
  y_out(:,flags) = y(2*N+1:end);
  wp_values(flags) = teichon_wp(y(1:N), y(N+1:2*N));
end

% Implementation:
while time<T;
  k = zeros([2*N+M 1]);
  stage_time = t;
  dt = t(step_number+1) - t(step_number);

  % RK stages
  for p = 1:length(rk.a);
    stage_time = time + dt*rk.c(p);

    f = epdiff_teichon_rhs(y(1:N),y(N+1:2*N), y(2*N+1:end));

    k = rk.a(p)*k + dt*f;

    y = y + rk.b(p)*k;

    if opt.wrap
      y(N+1:end) = wrap(y(N+1:end));
    end

  end

  % Update time, iteration count
  time = t(step_number+1);
  step_number = step_number+1;

  % Give output
  flags = abs(time-ts)<1e-12;
  if any(flags)
    a_out(:,flags) = y(1:N);
    b_out(:,flags) = y(N+1:2*N);
    y_out(:,flags) = y(2*N+1:end);
    
    wp_values(flags) = teichon_wp(y(1:N), y(N+1:2*N));
  end
end
