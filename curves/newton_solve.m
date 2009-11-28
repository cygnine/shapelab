function[x,converged] = newton_solve(z, x)
% newton_solve -- A Newton solver used by curvature_coefficients
%
% [x,converged] = newton_solve(z,x)
%
%     Solves the requisite Newton system as needed by
%     shapelab.curves.curvature_coefficients. 

persistent curv_func dcurv_func
if isempty(curv_func)
  from shapelab.curves import curvature_function as curv_func
  from shapelab.curves import curvature_function_jacobian as dcurv_func
end

N = length(z);
f = @(xx) curv_func(xx(1:(N-1)),xx(N:end),z);
df = @(xx) dcurv_func(xx(1:(N-1)), xx(N:end));

max_iter = 2000;
ftol = 1e-12;

N_iter = 0;

err = f(x);
while (norm(err)>ftol) & (N_iter < max_iter)
  dx = -inv(df(x))*err;
  x = x + dx;

  err = f(x);
  N_iter = N_iter + 1;
end

if norm(err)>ftol;
  fprintf('Failed to converge\n');
  converged = false;
else
  converged = true;
end
