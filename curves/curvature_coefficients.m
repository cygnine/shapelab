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
%          s = [1/2, 3], then z(2) = z(1) + int_0^{1/2} cos(kappa(s)) +
%                                                       i*sin(kappa(s)) ds
%
%                             z(3) = z(1) + int_0^{3}   cos(kappa(s)) +
%                                                       i*sin(kappa(s)) ds

persistent newton_solve eval_curvature_samples
if isempty(newton_solve)
  from shapelab.curves import eval_curvature_samples newton_solve
end

N = length(z);
% Translate all the points so z(1) is the origin
z = z(:) - z(1);

k = zeros([N-1 1]);
s = zeros([N-1 1]);

if N-2==0 % Line interpolation -- easy
  s = abs(z(2));
  k = atan(imag(z(2))/real(z(2)));
  return

else % Circle interpolation -- not too bad
  x = real(z(1:3));
  y = imag(z(1:3));
  a = abs(z(1:3));

  k(1) = (y(2)*a(3)^2 - y(3)*a(2)^2)/(x(2)*a(3)^2 - x(3)*a(2)^2);
  k(1) = atan2(y(2)*a(3)^2 - y(3)*a(2)^2, x(2)*a(3)^2 - x(3)*a(2)^2);

  k(2) = -2*(x(2)*sin(k(1)) - y(2)*cos(k(1)))/a(2)^2;

  s(1) = atan2(k(2)*x(2)+sin(k(1)), -k(2)*y(2) + cos(k(1)));
  s(1) = (s(1) - k(1))/k(2);

  s(2) = atan2(k(2)*x(3)+sin(k(1)), -k(2)*y(3) + cos(k(1)));
  s(2) = (s(2) - k(1))/k(2);
end

N = N-1;

% Otherwise: *sigh*
for q = 3:N  % Successively build up interpolations

  x = [k(1:q); s(1:q)]; % Use previous solution as guess
  x(end) = x(end-1);

  N_small_steps = 20;
  zend = eval_curvature_samples(x(1:(q-1)), x(end)+abs(z(q+1)-z(q)));
  z_ends = linspace(zend, z(q+1), N_small_steps+1);
  s_step = abs(z(q+1)-z(q))/N_small_steps;
  z_ends(1) = [];

  x(end) = x(end) + abs(z(q+1)-z(q));

  for qq = 1:N_small_steps
    ztemp = [z(1:q); z_ends(qq)];
    [x,converged] = newton_solve(ztemp, x);
    if not(converged)
      return
    end
  end

  k(1:q) = x(1:q);
  s(1:q) = x((q+1):end);

end
