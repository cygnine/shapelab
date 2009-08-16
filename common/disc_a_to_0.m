function[w] = disc_a_to_0(z,a)
% [w] = disc_a_to_0(z,a)
%
%     Uses a linear fractional map on the input z that takes the point a to 0
%     while preserving the unit circle.
%
%     The inverse of disc_0_to_a.

if abs(a)<1e-14
  w = z;
  return
elseif isinf(a)
  w = 1./z;
  return
end

global handles;
moebius = handles.shapelab.common.moebius;
H = abs(a)/a*[1 -a; -conj(a) 1];
w = moebius(z, H*abs(a)/a);
%w = abs(a)/a*(z-a)./(1-conj(a)*z);
