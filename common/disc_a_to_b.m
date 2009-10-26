function[w] = disc_a_to_b(z,a,b)
% [w] = disc_a_to_b(z,a,b)
%
%     Uses a linear fractional map taking a ---> b on the input z. This map
%     preserves the unit circle.
%
%     TODO: figure this out explicitly to ameliorate roundoff.

global packages;
da0 = packages.shapelab.common.disc_a_to_0;
d0b = packages.shapelab.common.disc_0_to_a;

w = d0b(da0(z,a),b);
