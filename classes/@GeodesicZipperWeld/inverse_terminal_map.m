function[z] = inverse_terminal_map(self, z, interior, slit_interior, slit_exterior, slit_interior_limbo, slit_exterior_limbo)
% inverse_terminal_map -- The inverse terminal map for zipper-type algorithms
%
% [z] = inverse_terminal_map(self, z, mapdata, [[interior, slit_interior,
% slit_exterior, slit_interior_limbo, slit_exterior_limbo]])
%
%     Performs the inverse of the termnal map for zipper-type algorithms. The
%     optional inputs are all indicial or boolean indexing arrays. 
%
%     `interior' -- any point not on the boundary of either the exterior or
%                   interior map
%     `slit_interior' -- any point on the boundary of the interior map
%     `slit_exterior' -- any point on the boundary of the exterior map
%     'slit_interior_limbo' -- point on slit_interior that need special
%       consideration due to complex square roots (they lie between z(N) and z(1)
%       on the original shape.
%     'slit_exterior_limbo' -- point on slit_exterior that need special
%       consideration due to complex square roots (they lie between z(N) and z(1)
%       on the original shape.
%
%     If boolean inputs are specified, this function does nothing to values of
%     z that are not indexed.


persistent csqrt
if isempty(csqrt)
  from shapelab.common import positive_angle_exponential as csqrt
end

sgn = -sign(self.winding_number);

if nargin<4
  z = csqrt(sgn*z, 1/2, 'cut_bias', sgn>0);
  z = self.moebius_maps.terminal_map.inv(z);
  %z = moebius_inverse(z, mapdata.moebius_maps.terminal_map);
else
  if nargin < 7
    slit_interior_limbo = false(size(slit_interior));
    slit_exterior_limbo = false(size(slit_exterior));
  end

  z(interior) = csqrt(sgn*z(interior), 1/2);
  z(slit_interior) = csqrt(sgn*z(slit_interior), 1/2, 'cut_bias', sgn>0);
  z(slit_exterior) = csqrt(sgn*z(slit_exterior), 1/2, 'cut_bias', sgn<0);

  z(slit_interior & ~slit_interior_limbo) = real(z(slit_interior & ~slit_interior_limbo));
  z(slit_exterior & ~slit_exterior_limbo) = real(z(slit_exterior & ~slit_exterior_limbo));

  allz = interior | slit_interior | slit_exterior;
  %z(allz) = moebius_inverse(z(allz), mapdata.moebius_maps.terminal_map);
  z(allz) = self.moebius_maps.terminal_map.inv(z(allz));
end
