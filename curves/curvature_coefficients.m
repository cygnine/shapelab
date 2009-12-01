function[k,s] = curvature_coefficients(z)
% curvature_coefficients -- returns a polynomial representation for curvature
%
% [k, s] = curvature_coefficients(z)
%
%     Interpolates the N complex-valued data points z with curve having
%     polynomial curvature with degree (N-2). The curvature is determined by
%     solving a 2*(N-2) component nonlinear system. The points z are assumed to
%     be ordered in the direction of interpolation. 
%   
%     The N-2 coefficients k are the monomial coefficient representation of the
%     primitive of the curvature. The N-1 coefficients s are the integral bounds
%     that interpolate to z.
%
%     E.g. k = [-1 2 3], then kappa(s) = -1 + 2*s + 3*s^2.
%          s = [0, 1/2, 3], then z(2) = z(1) + int_0^{1/2} cos(kappa(s)) +
%                                                       i*sin(kappa(s)) ds
%
%                             z(3) = z(1) + int_0^{3}   cos(kappa(s)) +
%                                                       i*sin(kappa(s)) ds

persistent newton_solve eval_curvature_samples circle_coefficients ...
           line_coefficients
if isempty(newton_solve)
  from shapelab.curves import eval_curvature_samples newton_solve ...
                              circle_coefficients line_coefficients
end

N = length(z);
% Translate all the points so z(1) is the origin
z = z(:) - z(1);

k = zeros([N-1 1]);
s = zeros([N-1 1]);

if N-2==0 % Line interpolation -- easy
  [k,s] = line_coefficients(z(1:2));
  return

else % Circle interpolation -- not too bad
  [k,s] = circle_coefficients(z(1:3));
end

N = N-1;

% Otherwise: *sigh*
for q = 3:N  % Successively build up interpolations
  [k,s] = extend_curve(k(1:q-1), s(1:q-1), z(1:q), z(q+1));
end
