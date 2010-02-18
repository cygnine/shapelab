function[H] = ba_driver(h, z)
% ba_driver -- Driver file for the Beurling-Ahlfors quasiconformal extension
%
% H = ba_driver(h, z)
%
%     The input h(.) will be extended from the real line to the complex points
%     z. This is effectively an implementation of (19) and (20) in [1]. No error
%     checking is performed to ensure that the input points z are in the upper
%     half-plane.
%
%      [1]: "Universal Teichmüller Space", F.P. Gardiner and W.J. Harvey, 2000
%           (arXiv print)

persistent gq spdiag r w N 
if isempty(gq)
  from speclab.orthopoly1d.jacobi.quad import gauss_quadrature as gq
  from labtools import spdiag

  N = 100;   % Number of quadrature points
  [r,w] = gq(N, 'alpha', 0, 'beta', 0);
  r = r.';
end

zsize = size(z);
z = z(:);

x = real(z);
y = imag(z);
Y = abs(y)/2;
Ym = spdiag(Y);

R1 = repmat(r, [length(z), 1]);
R2 = repmat(r, [length(z), 1]);

% Scale quadrature points to correct length
R1 = Ym*R1;
R2 = Ym*R2;
% Shift quadrature points to right centrices
R1 = R1 + repmat(x+Y, [1, N]);
R2 = R2 + repmat(x-Y, [1, N]);

% \int_{x}^{x+y} h(t) dt
H1 = Y.*(h(R1)*w);

% \int_{x-y}^{x} h(t) dt
H2 = Y.*(h(R2)*w);

H = 1./(2*y).*(H1 + H2);
H = H + i./y.*(H1 - H2);

H = reshape(H, zsize);
