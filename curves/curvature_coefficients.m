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

persistent newton_raphson curv_func dcurv_func
if isempty(newton_raphson)
  from labtools.rootfind import newton_raphson_system as newton_raphson
  from shapelab.curves import curvature_function as curv_func
  from shapelab.curves import curvature_function_jacobian as dcurv_func
end

N = length(z);
% Translate all the points so z(1) is the origin
z = z - z(1);

if N-2==0 % Line interpolation -- easy
  s = abs(z(2));
  k = atan(imag(z(2))/real(z(2)));
  return

%elseif N-2==1 % Circle interpolation -- not too bad
  x = real(z);
  y = imag(z);
  a = abs(z);

  k(1) = (y(2)*a(3)^2 - y(3)*a(2)^2)/(x(2)*a(3)^2 - x(3)*a(2)^2);
  k(1) = atan2(y(2)*a(3)^2 - y(3)*a(2)^2, x(2)*a(3)^2 - x(3)*a(2)^2);

  k(2) = -2*(x(2)*sin(k(1)) - y(2)*cos(k(1)))/a(2)^2;

  s(1) = atan2(k(2)*x(2)+sin(k(1)), -k(2)*y(2) + cos(k(1)));
  s(1) = (s(1) - k(1))/k(2);

  s(2) = atan2(k(2)*x(3)+sin(k(1)), -k(2)*y(3) + cos(k(1)));
  s(2) = (s(2) - k(1))/k(2);
  return
end

N = N-1;

% Otherwise: *sigh*
f = @(x) curv_func(x(1:N),x((N+1):end),z);
df = @(x) dcurv_func(x(1:N), x((N+1):end));

% Inital guess: linear connection from z(1) to z(end)
s = linspace(0, abs(z(end)-z(1)), N+1).';
s(1) = [];
k = zeros([N 1]);
k(1) = atan2(imag(z(end)-z(1)), real(z(end)-z(1)));

x = [k(:); s(:)];

% Newton's method on the above
% Rudimentary strategy:
max_iter = 1e3;
ftol = 1e-12;
baby_tol = 1e-1;

N_iter = 0;

err = f(x);
while (norm(err)>baby_tol) & (N_iter < max_iter)
  dx = -1/100*inv(df(x))*err;
  x = x + dx;

  err = f(x);
  N_iter = N_iter + 1;
end

while (norm(err)>ftol) & (N_iter < max_iter)
  dx = -inv(df(x))*err;
  x = x + dx;

  err = f(x);
  N_iter = N_iter + 1;
end

if norm(err)>ftol;
  x = NaN*x;
end

%x = newton_raphson(x0, f, df);

k = x(1:N);
s = x((N+1):end);
