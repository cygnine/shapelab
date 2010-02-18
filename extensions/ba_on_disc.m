function[H] = ba_on_disc(theta_int, theta_ext, w)
% ba_on_disc -- Evaluates the Beurling-Ahlfors extension for S^1 homeomorphisms
%
% H = ba_on_disc(theta_int, theta_ext, w)
%
%     Given point-values (theta_int, theta_ext) representing a
%     fingerprint on S^1, this function implements the BA extension on H of
%     the fingerprint extended (quasiperiodically) on R. It then maps these
%     values back to D.

persistent phi inv_phi ba
if isempty(phi)
  from shapelab.common import H_to_D_cover as phi
  from shapelab.common import D_to_H_cover as inv_phi
  from shapelab.extensions import ba_fingerprint as ba
end

H = inv_phi(ba(theta_int, theta_ext, phi(w)));
