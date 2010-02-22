function[mu] = beltrami_coefficient_driver(stuff, w)
% beltram_coefficient_driver -- Computes the Beltrami coefficient of the BA extension for a fingerprint
%
% mu = beltrami_coefficient_driver(stuff, w)
% 
%     Given a welding map as precomputed "stuff" (see
%     weld_primitive_setup), this computes the Beltrami coefficient
%     of the Beurling-Ahlfors extension at the locations z. For details about
%     the approximations used, see weld_primitive_evaluate and
%     weld_primitive_setup.

persistent inv_phi fprint_extend ba
if isempty(inv_phi)
  from shapelab.common import D_to_H_cover as inv_phi
  from shapelab.extensions import fingerprint_extend as fprint_extend
  from shapelab.extensions import ba_fingerprint_driver as ba
end

z = inv_phi(w);
x = real(z);
y = imag(z);

H = ba(stuff, z);
hx = fprint_extend(stuff,x); 
% 'h(x+y)'
hxpy = fprint_extend(stuff, x+y);
% 'h(x-y)'
hxmy = fprint_extend(stuff, x-y);

b = 3 + i;
bbar = conj(b);

% This is dh/dzbar divided by dh/dz:
mu = i*(b*hxpy + bbar*hxmy - 4*hx - 2*H)./...
     (b*hxpy - bbar*hxmy - 4*i*hx + 2*i*H);
%mu = i*(b*hxmy + bbar*hxpy - 4*hx - 2*conj(H))./...
%mu = conj(w)./(i*w).*mu;
mu = -mu.*w./conj(w);
