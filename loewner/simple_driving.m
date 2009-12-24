function[zc] = simple_driving(xi, T, z, varargin)
% simple_driving -- Low-accuracy solution to the Loewner equation
%
% zc = simple_driving(xi, T, {dt=0.01, visualize=false})
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

persistent input_schema
if isempty(input_schema)
  from labtools import input_schema
end

% rk data:
rka = [0 ...
    -567301805773/1357537059087 ...
    -2404267990393/2016746695238 ...
    -3550918686646/2091501179385 ...
    -1275806237668/842570457699];

rkb = [1432997174477/9575080441755 ...
    5161836677717/13612068292357 ...
    1720146321549/2090206949498 ...
    3134564353537/4481467310338 ...
    2277821191437/14882151754819];

rkc = [ 0.0 ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0 ...
    1.0];

inputs = {'dt', 'visualize'};
defaults = {0.01, false};

opt = input_schema(inputs, defaults, [], varargin{:});

N = ceil(T/opt.dt);
z = z(:);
zc = zeros(size(z));
niter = 0;

% Reflect xi so rk can work forward in time
xit = @(t) xi(T - t);
t = 0;

if opt.visualize
  plot_id = make_plots();
end

xin = xit(opt.dt);
% The first step is a backwards Euler step
%zc = backward_euler_step(z, 2*opt.dt, xin);
zc = z;

%zc = 1/2*(xin+z) + 1/2*sqrt((z + xin).^2 - 4*(xin.*z + opt.dt));
%flags = imag(zc)<0; % Sometimes need to reflect up from -H to H
%zc(flags) = conj(zc(flags));

realz = abs(imag(zc))<1e-13;

if opt.visualize
  update_plots()
end

while t < T
  ku = zeros(size(zc));
  kutemp = ku;
  stage_time = t;
  if t+opt.dt>T;
    opt.dt = T - t;
  end

  oldz = zc;
  flags = false(size(zc));
  % RK stages
  for p = 1:length(rka)
    stage_time = t + opt.dt*rkc(p);

    zc_rhs = -2./(zc - xit(stage_time));  % Need a negative sign: t is flipped around

    kutemp = ku;
    ku = rka(p)*ku + opt.dt*zc_rhs;

    newz = zc + rkb(p)*ku;

    % Do something different for evolutions that "jump" across the singularity
    newtime = t + opt.dt*rkc(p+1);
    flags = flags | sign(zc - xit(stage_time)).*sign(newz - xit(newtime)) < 0.5;
    flags = flags | not(isfinite(newz));

    zc = newz;
    realz = abs(imag(zc))<1e-13;
    %zc = zc + rkb(p)*ku;
  end

  flags = flags & realz;
  if any(flags)
    zc(flags) = backward_euler_step(oldz(flags), 2*opt.dt, xit(t));
  end

  realz = abs(imag(zc))<1e-13;

  if opt.visualize
    update_plots()
  end

  t = t+opt.dt;
  niter = niter+1;
end

function plot_id = make_plots()
  plot_id = plot(z, 'r.');
end

function [] = update_plots()
  set(plot_id, 'xdata', real(zc), 'ydata', imag(zc));
  axis([0, 2, 0, 2]);
  drawnow; pause(0.01)
end

end

function [x] = backward_euler_step(A,B,C)
% backward_euler_step -- Solves one backward Euler step
%
%     The equation is 
%
%          g_{n+1} = g_n + 2 dt/(g_n - xi(t_n)), rewritten as
%
%             A    = g_n + B /(g_n - C)

x = 1/2*(C+A + sqrt((C+A).^2 - 4*(B+A.*C)));
flags = imag(x)<0;
x(flags) = conj(x(flags));

end
