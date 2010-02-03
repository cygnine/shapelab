function[z] = inverse_terminal_map(z, mapdata, interior, slit_interior, slit_exterior)
% inverse_terminal_map -- The inverse terminal map for zipper-type algorithms
%
% [z] = inverse_terminal_map(z, mapdata, [[interior, slit_interior,
% slit_exterior]])
%
%     Performs the inverse of the termnal map for zipper-type algorithms. The
%     optional inputs are all indicial or boolean indexing arrays. 
%
%     `interior' -- any point not on the boundary of either the exterior or
%                   interior map
%     `slit_interior' -- any point on the boundary of the interior map
%     `slit_exterior' -- any point on the boundary of the exterior map
%
%     If boolean inputs are specified, this function does nothing to values of z
%     that are not indexed.


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
    z = moebius_inverse(z, mapdata.moebius_maps.terminal_map);
  else
    z(interior) = csqrt(sgn*z(interior), 1/2);
    z(slit_interior) = real(csqrt(sgn*z(slit_interior), 1/2, 'cut_bias', sgn>0));
    z(slit_exterior) = real(csqrt(sgn*z(slit_exterior), 1/2, 'cut_bias', sgn<0));
    allz = interior | slit_interior | slit_exterior;
    z(allz) = moebius_inverse(z(allz), mapdata.moebius_maps.terminal_map);
  end
otherwise
  error('Unrecognized algorithm type');
end
