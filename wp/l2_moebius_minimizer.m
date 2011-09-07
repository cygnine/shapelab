function[a,b,gradnorm] = l2_moebius_minimizer(theta, phi, varargin)
% l2_moebius_minimizer -- Computes the Moebius map connecting two data sets
%
% [a,b,gradnorm] = l2_moebius_minimzer(theta, phi, {cg_options})
%
%     Inputs theta and phi are sorted vectors of the same length. This function
%     seeks the two-parameter Moebius map defined by
%
%       phi = 2*atan(tan(theta/2)/(a*tan(theta/2) + b))
%
%     that minimizes the discrete L2 norm between phi and theta. Here phi and
%     theta are angles over [0, 2*pi). The map has only two parameters (instead
%     of three) because theta=0 is assumed to map to phi=0. 
%
%     Any optional inputs are passed along to labtools.minimization.nonlinear_cg
%     -- see that function for explanation of optional inputs.

persistent mymod nonlinear_cg
if isempty(mymod)
  from labtools import interval_wrap as mymod
  from labtools.minimization import nonlinear_cg
end

moebius_eval = @(a,b,t) mymod(2*atan(tan(t/2)./(a*tan(t/2)+b)), [0, 2*pi]);

N = length(theta);
if length(phi) ~= N
  error('Inputs theta and phi must have the same length');
end
if any(diff(theta)<0) | any(diff(phi)<0)
  error('Inputs theta and phi must be sorted');
end

L_diff = @(a,b,t) moebius_eval(a,b,t) - phi;
L_eval = @(a,b,t) sum(abs(L_diff(a, b, t)).^2);
% The actual objective function of one variable:
L_obj = @(z) L_eval(z(1), z(2), theta);

phi_grad_a = @(a,b,t) -2./(1 + (a+b*cot(t/2)).^2);
phi_grad_b = @(a,b,t) -2./((a + b*cot(t/2)).*(a*tan(t/2) + b) + tan(t/2));
L_grad = @(a,b,t) [2*L_diff(a,b,t).'*phi_grad_a(a,b,t); ...
                       2*L_diff(a,b,t).'*phi_grad_b(a,b,t)];
% The actual gradient of the objective (function of one variable):
L_objgrad = @(z) L_grad(z(1), z(2), theta);

% The hard work:
[z, gradnorm] = nonlinear_cg(L_obj, L_objgrad, [0; 1], varargin{:});
a = z(1);
b = z(2);
