function[w] = disc_0_to_a(z,a)
% [w] = disc_0_to_a(z,a)
%
%     Uses a linear fractional map on the input z that takes the point 0 to a
%     while preserving the unit circle.
%
%     The inverse of disc_a_to_0.

if abs(a)<1e-14
  w = z;
  return
elseif isinf(a)
  w = 1./z;
  return
end

b = abs(a);
w = a/b*(z+b)./(1+b*z);
