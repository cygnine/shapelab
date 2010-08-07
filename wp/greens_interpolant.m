function[f] = greens_interpolant(theta, v, z, suppress_kernel)
% greens_interpolant -- Interpolates samples using the Green's function
%
% f = greens_interpolant(theta, v, z,{{suppress_kernel=true}})
%
%     Given samples (theta,v) of a velocity field v, this function constructs
%     the Green's function interpolant and evaluates it at the locations z \in
%     S^1. The locations theta must be distinct (no checking is performed).
%     The optional flag suppress_kernel (default true) tells the function
%     whether or not to push elements in the kernel through the inversion of G.
%
%     If suppress_kernel is true, the steps for the procedure are the following:
%        1.) Compute the interpolatory Vandermonde-like matrix G. Note
%            that since L is positive, G is positive as well. Direct eig solver
%            is used, but since G is spd, mayhaps iterative methods are better.
%        2.) Compute the eig-decomposition of G. Since dim(ker(L)) = 3 and is
%            given by (1, sin, cos), then we effectively orthogonalize v against
%            the kernel and then invert G to obtain coefficients.
%        3.) Components of the kernel are propagated unchanged and added to the
%            interpolant.
%        4.) Compute the resulting interpolant at z.
%
%     When suppress_kernel is set to false, nothing fancy is done; G is inverted
%     using a Cholesky decomposition.

persistent greens_function L_kernel_vandermonde greens_coefficients
if isempty(greens_function)
  from shapelab.wp import greens_function greens_coefficients
  L_kernel_vandermonde = @(x) [ones([length(x) 1]) sin(x(:)) cos(x(:))];
end

if nargin<4
  suppress_kernel = 2;
end

[u,u_kernel] = greens_coefficients(theta, v, suppress_kernel);

%% Form G
%[thetai, thetaj] = meshgrid(theta, theta);
%G = greens_function(thetai-thetaj);
%
%if suppress_kernel
%
%  [Gv,Gd] = eig(G);  % G is spd, this should be fast, stable 
%
%  % Orthogonalize v against 1, cos, sin.
%  % First compute G-spectral representation of (1, cos, sin):
%  L_kernel_spectrum = Gv'*L_kernel_vandermonde(theta);
%
%  % Now:
%  % 1.) kernel is orth to non-kernel in space
%  % 2.) Gv is unitary ====> kernel spectrum is orth to non-kernel spectrum
%  % 3.) Ergo, orthogonalize v-spectrum to kernel spectrum
%
%  % Must orthogonalize vectors in kernel spectrum first, form thin QR
%  [q,r] = qr(L_kernel_spectrum,0);
%
%  v_spectrum = Gv'*v(:);      % G-spectrum of v
%  u_kernel = q'*v_spectrum;  % non-normalized coeffs of (1,sin,cos) in v
%  v_kernel_spectrum = q*u_kernel;  % G-spectrum of u_kernel
%
%  % Now just invert G as usual
%  u = Gv*inv(Gd)*(v_spectrum - v_kernel_spectrum);  % non-kernel part
%else
%
%  % cholesky is probably better (quicker, more stable) if G is spd
%  R = chol(G);
%  u = R\(R'\v);
%end

% Compute interpolant evaluations for elements orthogonal to the kernel
[thetai, thetaj] = meshgrid(theta, z);
f = greens_function(thetai-thetaj)*u;  % non-kernel part

if suppress_kernel>0
  % Now add in residual parts in the kernel
  f = f + L_kernel_vandermonde(z)*u_kernel;
end

f = reshape(f, size(z));
