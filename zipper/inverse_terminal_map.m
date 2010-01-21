function[z] = inverse_terminal_map(z, mapdata, interior, slit_interior, slit_exterior)
% inverse_terminal_map -- The inverse terminal map for zipper-type algorithms
%
% [z] = inverse_terminal_map(z, mapdata)

persistent moebius_inverse csqrt
if isempty(moebius_inverse)
  from shapelab.common import moebius_inverse
  from shapelab.common import positive_angle_exponential as csqrt
end

switch lower(mapdata.type)
case {'geodesic', 'slit', 'zipper_weld', 'loewner'}
  sgn = -sign(mapdata.winding_number);
  if nargin<3
    z = csqrt(sgn*z, 1/2, 'cut_bias', sgn>0);
  else
    z(interior) = csqrt(sgn*z(interior), 1/2);
    z(slit_interior) = csqrt(sgn*z(slit_interior), 1/2, 'cut_bias', sgn>0);
    z(slit_exterior) = csqrt(sgn*z(slit_exterior), 1/2, 'cut_bias', sgn<0);
  end
  z = moebius_inverse(z, mapdata.moebius_maps.terminal_map);
otherwise
  error('Unrecognized algorithm type');
end
