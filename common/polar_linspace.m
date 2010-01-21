function[z] = polar_linspace(N,M,varargin)
% [z] = polar_linspace(N,M,{r0=0, r1=1, theta0=0, theta1=2*pi})
%
%     Creates a 2-D 'linspace' of N*M data points. Creates a polar coordinate
%     grid with M equispaced points in the theta direction between theta0 and
%     theta1, and N equispaced points in the r direction between r0 and r1.
%
%     Unlike Matlab's lin/logspace functions, this one returns a column vector.
%
%     If r0=0, this function isn't smart enough to see that and puts M points at
%     0.

persistent input_schema
if isempty(input_schema)
  from labtools import input_schema
end

inputs = {'r0', 'r1', 'theta0', 'theta1'};
defaults = {0,1,0,2*pi};
opt = input_schema(inputs, defaults, [], varargin{:});

z = zeros([N*M 1]);

r = linspace(opt.r0, opt.r1, N).';
if abs(mod(opt.theta1,2*pi) - opt.theta0)<1e-12
  theta = linspace(opt.theta0, opt.theta1, M+1).';
  theta = theta(1:M);
else
  theta = linspace(opt.theta0, opt.theta1, M).';
end
theta = exp(i*theta);

for q = 1:M
  i1 = 1 + (q-1)*N;
  i2 = q*N;
  z(i1:i2) = r*theta(q);
end

z = reshape(z, [N M]);
