function[f] = epdiff_teichon_rhs(a, b, y)
% epdiff_teichon_rhs -- evaluates the RHS of epdiff
%
% f = epdiff_teichon_rhs(a, b, y)
%
%     Evaluates the right-hand side of the Euler-Poincare differential equation
%     on Diff(S^1) for a Teichon ansatz. The vectors a and b have length N,
%     which corresponds to N teichons with strengths a, and positions b. The
%     evolution of these parameters is given by a differential equation; this
%     function evaluates the right-hand side of this differential equation. The
%     system of ODEs is 
%
%      a_k' = -a_k \sum_{j=1}^N G'(b_k - b_j) a_j
%
%      b_k' = \sum_{j=1}^N a_j G(b_k - b_j)
%
%      y_k' = \sum_{j=1}^N a_j G(y_k - b_j)

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
%Gp1 = greens_function_derivative(temp1-temp2);
Gp1 = greens_function(temp,1);

[temp1,temp2] = ndgrid(y, b);
G2 = greens_function(temp1-temp2);

f = zeros([2*N+M 1]);
f(1:N) = -spdiag(a)*(Gp1*a);
f((N+1):2*N) = G1*a;
f(2*N+1:end) = G2*a;
