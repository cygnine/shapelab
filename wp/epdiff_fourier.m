function[v_out, wp_values] = epdiff_fourier(v0, m0, ts, varargin)
% epdiff_fourier - Solves epdiff on Diff(S^1) using a Fourier method
%
% v = epdiff_fourier(v0, ts, {interpolation=cubic, filter='', N=256,dt=0.01})
%
%     Solves epdiff on Diff(S^1) using a global Fourier approximation.  This
%     differential equation is given by 
%
%      m_t + (m v)_\theta + v_\theta m = 0,
%      m = L v = - H * (v_{3\theta} + v_\theta),   \theta \in [0, 2*pi)
%
%     and H is the periodic Hilbert transform.  v0 and m0 are function handles
%     of one variable for evaluating the theta dependence of the initial data.
%     ts is a vector of times for which the output velocity is requested. The
%     spatial method is Galerkin.

persistent rk strict_inputs ffft iffft fft_setup fgq integer_range
persistent L Linv multiply wp_norm fourier_d
if isempty(rk)
  from odesolve.coeffs import lserk4
  from labtools import strict_inputs
  from speclab.fourier.fft import ffft_online as ffft
  from speclab.fourier.fft import iffft_online as iffft
  from speclab.fourier.fft import ffft_overhead as fft_setup
  from speclab.fourier import gauss_quadrature as fgq
  from speclab.fourier import coefficient_derivative_function as fourier_d
  from speclab.common import integer_range

  from speclab.fourier import coefficient_convolution_function as multiply
  from shapelab.diffs1 import L_fourier_coeffs as L
  from shapelab.diffs1 import invert_L_fourier_coeffs as Linv
  from shapelab.diffs1 import wp_norm_fourier_coeffs as wp_norm

  rk = lserk4();
end

inputs = {'filter', 'N', 'dt'};
defaults = {[], 256, 0.01};
opt = strict_inputs(inputs, defaults, [], varargin{:});
if isempty(opt.filter)
  opt.filter = ones([opt.N 1]);
end

fft_overhead = fft_setup(opt.N);

% First construct temporal vector
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

v_out = zeros([opt.N, length(ts)]);
[theta,w] = fgq(opt.N, 'shift', pi);
quotient_space = find(abs(integer_range(opt.N))<=1);
dtheta = fourier_d(opt.N);
fmultiply = multiply(opt.N);

% Get initial m, v (in spectral space)
v = ffft(v0(theta), fft_overhead);

if isempty(m0)
  m = L(v);
else
  % Mayhaps these aren't necessary?
  m = ffft(m0(theta), fft_overhead);
  assert(norm(Linv(m) - v)<1e-6, 'Error: initial data doesn''t satisfy m = L v');
end

time = 0;
step_number = 1;
wp_values = zeros(size(ts));

flags = abs(time-ts)<1e-12;
if any(flags)
  v_out(:,flags) = v;
  wp_values(flags) = wp_norm(v);
end

% Implementation:
while time<T;
  ku = zeros(size(m));
  stage_time = t;
  dt = t(step_number+1) - t(step_number);

  % RK stages
  for p = 1:length(rk.a);
    stage_time = time + dt*rk.c(p);

    u_rhs = -fmultiply(dtheta(v),m) - dtheta(fmultiply(v, m));

    ku = rk.a(p)*ku + dt*u_rhs;
    m = m + rk.b(p)*ku;
    % Project out any components in -1, 0, 1:
    m(quotient_space) = 0;
    v = Linv(m);
  end

  % Update time, iteration count
  time = t(step_number+1);
  step_number = step_number+1;

  % Give output
  flags = abs(time-ts)<1e-12;
  if any(flags)
    v_out(:,flags) = v;
    wp_values(flags) = wp_norm(v);
  end
end

v_out = iffft(v_out, fft_overhead);
