function[phi] = theta_automorphism(theta, H)
% theta_automorphism -- An automorphism of [0, 2*pi] from PSL_2(R)
%
% phi = theta_automorphism(theta,H)
%
%     Given H \in PSL_2(R), this function effectively implements the projective
%     transformation of the real line as dictated by H, except that it
%     identifies points on \theta \in [0, 2*pi) as points on the real line via
%     the mapping x = tan(theta/2).

persistent moebius wrap
if isempty(wrap)
  from shapelab.common import moebius
  from labtools import interval_wrap as wrap
end

phi = tan(theta/2);
phi = moebius(phi,H);
phi = wrap(2*atan(phi), [0, 2*pi]);
