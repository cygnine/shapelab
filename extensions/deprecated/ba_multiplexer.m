function[H] = ba_multiplexer(h, z)
% ba_multiplexer -- `Multiplexer' for the Beurling-Ahlfors quasiconformal extension
%
% H = ba_multiplexer(h, z)
%
%     The input h(.) will be extended from the real line to the complex points
%     z. This file does not perform any considerable computation -- it simply
%     manages calls to different functions based on the inputs z.

persistent extend
if isempty(extend)
  from shapelab.extensions import ba_driver as extend
end

zsize = size(z);
iz = imag(z);

% The 3 regions: real line, H, and lower half-plane
reals = (iz==0);
Hs = (iz>0);
nHs = not(reals | Hs);

H = zeros(size(z));

% Different things for each region
if any(reals)
  H(reals) = z(reals);
end
if any(Hs)
  H(Hs) = extend(h, z(Hs));
end
if any(nHs)
  H(nHs) = conj(extend(h, conj(z(nHs))));
end

H = reshape(H, zsize);
