function[J,eprhs] = epdiff_teichon_variation_rhs(a, b, y)
% epdiff_teichon_variation_rhs -- evaluates the Jacobian of the RHS of epdiff
%
% [J, eprhs] = epdiff_teichon_variation_rhs(a, b, y)
%
%     Evaluates the Jacobian right-hand side of the Euler-Poincare differential
%     equation on Diff(S^1) for a Teichon ansatz. The vectors a and b have
%     length N, which corresponds to N teichons with strengths a, and positions
%     b. The vector y is a vector of locations that are pushed by the
%     teichon-induced velocity field.  The evolution of these parameters is
%     given by a differential equation; this function evaluates the Jacobian of
%     the right-hand side of this differential equation with respect to a, b,
%     and y.
%
%     Since the rhs of the regular epdiff is already computed (eprhs), might as
%     well return that as well.

persistent greens_function spdiag
if isempty(greens_function)
  from shapelab.wp import greens_function 
  from labtools import spdiag
end

N = length(a);
M = length(y);
%[temp1, temp2] = meshgrid(b, b);
[temp1, temp2] = ndgrid(b, b);
temp = temp1-temp2;
G1 = greens_function(temp);
Gp1 = greens_function(temp, 1);
Gpp1 = greens_function(temp, 2);
%Gp1 = greens_function_derivative(temp1-temp2);
%Gpp1 = greens_function_derivative(temp1-temp2, 2);
A = spdiag(a);

[temp1,temp2] = ndgrid(y, b);
temp = temp1 - temp2;
G2 = greens_function(temp);
Gp2 = greens_function(temp, 1);

eprhs = zeros([2*N+M 1]);
eprhs(1:N) = -A*(Gp1*a);
eprhs((N+1):2*N) = G1*a;
eprhs(2*N+1:end) = G2*a;

J = zeros(2*N+M);
J(1:N,1:N) = -spdiag(Gp1*a) - A*Gp1;
J(1:N,(N+1):2*N) = A*Gpp1*A - spdiag(A*Gpp1*a);
J((N+1):2*N,1:N) = G1;
J((N+1):2*N,(N+1):2*N) = -Gp1*A + spdiag(Gp1*a);

% a and b aren't influeced by positions y:
% J(1:2*N,2*N+1:end) = 0;

% y is influenced by a, b, y:
J(2*N+1:end,1:N) = G2;
J(2*N+1:end,N+1:2*N) = -Gp2*A;
J(2*N+1:end,2*N+1:end) = diag(Gp2*a);
