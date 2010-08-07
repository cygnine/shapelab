function[u, u_kernel] = greens_coefficients(theta, v, suppress_kernel)
% greens_coefficients -- Computes Green's representation coefficients
%
% f = greens_coefficients(theta, v, {{supporess_kernel=1}})
%
%     Given samples (theta,v) of a velocity field v, this function constructs
%     the Green's function coefficients u such that the reconstructed function
%     \sum_j u(j) G(theta-theta(j)) is an interpolant of v at theta(j).
%
%     By default, the kernel of the operator L is `suppressed'; this behavior
%     can be overwritten by specifying the optional flag `suppress_kernel' to
%     e.g. 0. Possibilities are given below:
%
%     suppress_kernel = 0 : The inner product matrix is directly inverted.
%     suppress_kernel = 1 : Any portion of the input that `appears' to be in the
%                           kernel of the operator is subtracted out before
%                           inversion.
%     suppress_kernel = 2 : A portion of the kernel is subtracted out to
%                           minimize the norm.

persistent greens_function L_kernel_vandermonde
if isempty(greens_function)
  from shapelab.wp import greens_function
  L_kernel_vandermonde = @(x) [ones([length(x) 1]) sin(x(:)) cos(x(:))];
end

if nargin<3
  suppress_kernel = 2;
end

% Form G
[thetai, thetaj] = meshgrid(theta, theta);
G = greens_function(thetai-thetaj);

if suppress_kernel==1

  [Gv,Gd] = eig(G);  % G is spd, this should be fast, stable 

  % Orthogonalize v against 1, cos, sin.
  % First compute G-spectral representation of (1, cos, sin):
  L_kernel_spectrum = Gv'*L_kernel_vandermonde(theta);

  % Now:
  % 1.) kernel is orth to non-kernel in space
  % 2.) Gv is unitary ====> kernel spectrum is orth to non-kernel spectrum
  % 3.) Ergo, orthogonalize v-spectrum to kernel spectrum

  % Must orthogonalize vectors in kernel spectrum first, form thin QR
  [q,r] = qr(L_kernel_spectrum,0);

  v_spectrum = Gv'*v(:);      % G-spectrum of v
  u_kernel = q'*v_spectrum;  % non-normalized coeffs of (1,sin,cos) in v
  v_kernel_spectrum = q*u_kernel;  % G-spectrum of u_kernel

  % Now just invert G as usual
  u = Gv*inv(Gd)*(v_spectrum - v_kernel_spectrum);  % non-kernel part
  u_kernel = r\u_kernel;
elseif suppress_kernel==2
  R = chol(G);

  L_kernel_vand = L_kernel_vandermonde(theta);
  L_kernel_rep = R'\L_kernel_vand; % inv(R')*L_kernel
  v_rep = R'\v(:); % inv(R')*v;

  % The following two lines are the result of a straightforward minimization
  % procedure to determine how much of the kernel to subtract out.
  mat = L_kernel_rep'*L_kernel_rep;
  v_kernel_rep = L_kernel_rep'*v_rep;

  u_kernel = mat\v_kernel_rep;  % coefficients of (1, sin, cos) to subtract out
  u = R\(v_rep - L_kernel_rep*u_kernel);
else

  % cholesky is probably better (quicker, more stable) if G is spd
  R = chol(G);
  u = R\(R'\v(:));
  u_kernel = [];
end
