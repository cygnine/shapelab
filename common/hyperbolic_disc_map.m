function[w] = hyperbolic_disc_map(z,a)
% hyperbolic_disc_map -- a disc-preserving moebius map on hyperbolic coordinates
%
% w = hyperbolic_disc_map(z,a)
%
%     Implements the disc-preserving Moebius transform taking a ---> 0. However,
%     the input z = rho*exp(i*theta), where theta is the usual angle but rho is
%     the hyperbolic distance from 0. Returns w in the same format: w =
%     omega*exp(i*phi), where omega is the hyperbolic distance from 0.

persistent da0 he
if isempty(da0)
  from shapelab.common import disc_a_to_0 as da0
  from shapelab.common import hyperbolic_to_euclidean as he
end

rho = abs(z);
ze = he(z);
etheta = z./rho;
etheta(isnan(etheta)) = 0;

%eprho = exp(rho);
emrho = exp(-2*rho);

temp1 = abs((1-conj(a)*etheta) + emrho.*(1+conj(a)*etheta));
temp2 = abs((etheta-a) - emrho.*(etheta+a));

omega = 1/2*log(temp1 + temp2)-1/2*log(temp1-temp2);

% Just use the Euclidean map for the angle
ephi = exp(i*angle(da0(ze,a)));

w = omega.*ephi;
