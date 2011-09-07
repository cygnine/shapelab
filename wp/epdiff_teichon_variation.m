function[a, b, y, J, wp] = epdiff_teichon_variation(a0, b0, y0,T, varargin)
% epdiff_teichon_variation -- Evolves epdiff using a teichon ansatz
%
% [a, b, y, J, wp] = epdiff_teichon_variation(a0, b0, T, {dt=0.01,
%                                                     wrap=false})
%
%     Solves epdiff on elements in Diff(S^1) using a teichon ansatz.  The
%     teichons have velocity representation:
%
%      v(theta) = sum_j a_j G(theta-b_j) + sum_q k_q phi_q(theta),
%
%     Inputs a0 and b0 are the initial strengths and locations of the teichons.
%
%     Evolution proceeds using an ODE for the positions and locations of the
%     teichons; a 5-stage 4th order SSP LSERK method is used. The outputs are
%     returned at the terminal time T. The initial time is assumed to be t=0.
%
%     where G is the WP-Green's function, and phi_q for q = [1,2,3] is a basis
%     for the kernel of the WP norm (see shapelab.wp.wp_minimal_metric and
%     shapelab.wp.kernel_basis). Note that whereas a_j and b_j depend on time,
%     k_q does not.
%
%     The output J is the (2N+M x 2N+M) Jacobian of [a,b,y]' with respect to the
%     initial data [a0, b0, y0].
%
%     If the optional input wrap is set to true, then the locations b of the
%     teichons and positions y are wrapped to the [0, 2*pi) interval.

persistent strict_inputs epdiff_rhs
persistent rk interval_wrap wrap
persistent greens_function
if isempty(strict_inputs)
  from labtools import strict_inputs interval_wrap
  from shapelab.wp import epdiff_teichon_variation_rhs as epdiff_rhs
  from shapelab.wp import greens_function

  from odesolve.coeffs import lserk4
  rk = lserk4();
  wrap = @(t) interval_wrap(t, [0, 2*pi]);
end

N = length(a0);
M = length(y0);
if N < 4
  fprintf('With less than four samples, the evolution is trivial...aborting\n');
  a = [];
  b = [];
  y = [];
  J = [];
  wp = [];
end

opt = strict_inputs({'dt', 'wrap'}, {0.01, false}, [], varargin{:});

b0 = b0(:); a0 = a0(:);

% Randomly: dt = 0.01
dt = opt.dt;
t = 0:dt:T;
if abs(T-t(end))>1e-13
  t(end+1) = T;
end

% Initial data
b = b0;
a = a0;
y = y0;
J = eye(2*N+M);
if opt.wrap
  b = wrap(b);
end

time = 0;
step_number = 1;
wp = 0;

% Implementation:
while time<T;
  ka = zeros(size(a));
  kb = zeros(size(b));
  ky = zeros(size(y));
  kJ = zeros(size(J));
  dt = t(step_number+1) - t(step_number);

  % RK stages
  for p = 1:length(rk.a);
    stage_time = time + dt*rk.c(p);

    [Jrhs, eprhs] = epdiff_rhs(a, b, y);

    ka = rk.a(p)*ka + dt*eprhs(1:N);
    kb = rk.a(p)*kb + dt*eprhs((N+1):2*N);
    ky = rk.a(p)*ky + dt*eprhs(2*N+1:end);
    kJ = rk.a(p)*kJ + dt*Jrhs*J;

    a = a + rk.b(p)*ka;
    b = b + rk.b(p)*kb;
    y = y + rk.b(p)*ky;
    J = J + rk.b(p)*kJ;

    if opt.wrap
      b = wrap(b);
      y = wrap(y);
    end
  end

  % Update time, iteration count
  time = time + dt;
  step_number = step_number+1;

end
[temp1, temp2] = meshgrid(b,b);
G = greens_function(temp1-temp2);
wp = sqrt(a'*G*a);
