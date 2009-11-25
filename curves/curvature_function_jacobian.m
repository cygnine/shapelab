function[J] = curvature_function_jacobian(k,s)
% curvature_function_jacobian -- evaluates nonlinear function jacobian
%
% J = curvature_function_jacobian(k,s)
%
% If length(s) == length(k), then this function evaluates a square matrix of
% size 2*length(s). The matrix is the Jacobian of the function defined by
% curvature_function with respect to the variables k and s.

persistent trig_integral2
if isempty(trig_integral2)
  from shapelab.curves import trig_integral2
end

N = length(s);
J = zeros(2*N);

cells = zeros([N 2]);
cells(:,1) = [0; s(1:N-1)];
cells(:,2) = s(:);
%%% Now construct Jacobian

% First N columns: derivatives with respect to k.
[cint,sint] = trig_integral2(k,cells);

temp = [-sint; cint];
re_index = reshape(1:(2*N), [N 2]).';
re_index = re_index(:);
J(:,1:N) = temp(re_index,:);

% Last N columns: derivaties with respect to s
% Easier to compute, harder to code
s = s(:).';
k = flipud(k(:));
p = polyval(k,s);

cint = cos(p);
sint = sin(p);
submat = zeros([2*N N]);

for q = 1:N
  r1 = 2*q-1;
  r2 = 2*q+2;

  if q==N
    submat(r1:end,q) = [cint(end); sint(end)];
  else
    submat(r1:r2,q) = [cint(q); sint(q); -cint(q); sint(q)];
  end
end

J(:,(N+1):end) = submat;
