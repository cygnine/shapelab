function[w] = disc_a_to_b(z,a,b)
% [w] = disc_a_to_b(z,a,b)
%
%     Uses a linear fractional map taking a ---> b on the input z. This map
%     preserves the unit circle.
%
%     TODO: figure this out explicitly to ameliorate roundoff.

persistent da0 d0b
if isempty(da0)
  from shapelab.common import disc_a_to_0 as da0
  from shapelab.common import disc_0_to_a as d0b
end

w = d0b(da0(z,a),b);
