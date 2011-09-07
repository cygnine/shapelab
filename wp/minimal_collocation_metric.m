function[H, K, Gi, P] = minimal_metric(theta)
% minimal_metric -- Computes the metric matrix for the `minimal' WP lift
%
% H = minimal_metric(theta)
% [H, K] = minimal_metric(theta)
% [H, K, Ginv, P] = minimal_metric(theta)
%
%     Given a vector of locations theta on S^1 = [0, 2*pi), this function
%     computes the (non-invertible) metric matrix H such that if v is a vector
%     of evaluations at theta, v.'*H*v returns the WP-norm (squared) of the
%     unique H^{3/2} function that (a) interpolates the data v, and (b) has
%     minimal WP norm.
%
%     The precise function from H^{3/2} that is this minimizer is given by 
%
%      w = \sum_{m=1}^M (H*v)(m) * G(. - theta(m)) + 
%          \sum_{m=1}^3 (K*v)(m) * phi(m),
%
%     where phi = [1, sin, cos], a basis for the kernel of WP norm, and G(theta)
%     is the Green's function under the WP metric.
%
%     The matrix H is given by Ginv*P, where P is a (Ginv-orthogonal) projection
%     matrix and Ginv is the inverse of the Gram matrix for Green's functions
%     centered at theta in the WP inner product.

persistent greens_function kernel_basis
if isempty(greens_function)
  from shapelab.wp import kernel_basis
  from shapelab.wp.teichon import greens_matrix
end

N = length(theta);
G = greens_matrix(theta);

U = kernel_basis(theta);

R = chol(G);
Ri = inv(R);
Gi = Ri*Ri';

K = inv(U'*Gi*U)*U'*Gi; % K = "kernel coefficients"

P = eye(N) - U*K;
H = Gi*P;
