function[dphi] = theta_automorphism_derivative(theta, H)
% theta_automorphism_derivative -- Derivative of theta_automorphism
%
% dphi = theta_automorphism_derivative(theta,H)
%
%     This computes the derivative of the isomorphic map theta_automorhpism.

persistent moebius moebius_derivative wrap
if isempty(wrap)
  from shapelab.common import moebius moebius_derivative
  from labtools import interval_wrap as wrap
end

phi = tan(theta/2);
phi = moebius(phi,H);
phi = wrap(2*atan(phi), [0, 2*pi]);

x = tan(theta/2);
y = moebius(x, H);

dphi = 2./(1+y.^2).* ...
       moebius_derivative(x, H).* ...
       1/2.*sec(theta/2).^2;
