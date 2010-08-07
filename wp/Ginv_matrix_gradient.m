function[Ginv, factors, Vs] = Ginv_matrix_gradient(theta)
% Ginv_matrix_gradient -- gradient of the Green's function matrix
%
% [Ginv, factors, Vs] = Ginv_matrix_gradient(theta)
%
%     The Green's function matrix G is a symmetric positive definite matrix,
%     whose (i,j) entry is the evaluation Gf(theta(i) - theta(j)), where Gf is
%     the Green's function. This function returns the quantities necessary to
%     construct the rank-2 derivative of inv(G) with respect to theta(j) for
%     all j. 
%
%     d G/d theta(j) = Vs{j}*diag([factors(j), -factors(j)])*Vs{j}'
%
%     hence factors(j) is the magnitude of the eigenvalue of dG/dtheta(j), and
%     Vs{j} contain the eigenvectors.

persistent Gf Gfd
if isempty(Gf)
  from shapelab.wp import greens_function as Gf
  from shapelab.wp import greens_function_derivative as Gfd
end

N = length(theta);
[thetai, thetaj] = meshgrid(theta);
thetas = thetaj - thetai;
G = Gf(thetas);
Ginv = inv(G);

G_d_theta = Gfd(thetas);
factors = sqrt(sum(G_d_theta.^2, 1)).';

Vs = cell([N 1]);

for n = 1:N
  z = G_d_theta(:,n)/(sqrt(2)*factors(n));
  en = zeros([N 1]); en(n) = 1;

  Vs{n} = Ginv*[z+en/sqrt(2), -z+en/sqrt(2)];
end
