function[w] = disc_0_to_a(z,a)
% [w] = disc_0_to_a(z,a)
%
%     Uses a linear fractional map on the input z that takes the point 0 to a
%     while preserving the unit circle.
%
%     The inverse of disc_a_to_0.

persistent moebius
if isempty(moebius)
  from shapelab.common import moebius
end

if abs(a)<1e-14
  w = z;
  return
elseif isinf(a)
  w = 1./z;
  return
end

%b = abs(a);
%H = [1 b; b 1];
H = abs(a)/a*[1 -a; -conj(a) 1];
H = [H(2,2) -H(1,2); -H(2,1) H(1,1)];
w = moebius(z,H);
%w = moebius(z, H*a/b);

%b = abs(a);
%w = a/b*(z+b)./(1+b*z);
