function[Hi, K] = wp_minimal_metric_inverse(theta)
% wp_minimal_metric -- Computes the inverse metric matrix for the `minimal' WP extension
%
% Hi = wp_minimal_metric_inverse(theta)
% [Hi, K] = wp_minimal_metric_inverse(theta)
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

persistent greens_function kernel_basis
if isempty(greens_function)
  from shapelab.wp import greens_function kernel_basis
  %L_kernel_vandermonde = @(x) [ones([length(x) 1]) sin(x(:)) cos(x(:))];
end

N = length(theta);
[thetai, thetaj] = meshgrid(theta, theta);
G = greens_function(thetai-thetaj);
U = kernel_basis(theta);

R = chol(G);
Ri = inv(R);
Gi = Ri*Ri';

K = inv(U'*Gi*U)*U'*Gi; % K = "kernel coefficients"

H = eye(N) - U*K;
%H = Gi*H;
Hi = H*G;
