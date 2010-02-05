function[vertices_int, vertices_ext, z, assumptions] = force_inverse_terminal_map(vertices_int, vertices_ext, z, assumptions, interior, slit_interior, slit_exterior)
% force_inverse_terminal_map -- The inverse terminal map for zipper-type algorithms
%
% [vertices_int, vertices_ext, z, assumptions] = force_inverse_terminal_map(vertices_int, vertices_ext, z, assumptions, [[interior, slit_interior, slit_exterior]])
%     Performs the necessary operations so that an inverse terminal map for
%     zipper-type algorithms is implemented. The only data available are the
%     interior and exterior vertices, quantities from the welding map. This is a
%     canonical inverse terminal map, meaning that it is assumed that:
%
%       vertices_int(end) == vertices_ext(end) == 0,
%       vertices_int(1) == vertices_ext(1) == Inf,
%
%     The terminal map operation is performed for every element in z. See
%     inverse_terminal_map for explanations of the boolean indexing vectors
%     `interior', `slit_interior', and `slit_exterior'.
%
%     The mandatory input `assumptions' is a struct containing the
%     assumptions about the map.

persistent moebius_inverse csqrt find_moebius invert_moebius
if isempty(moebius_inverse)
  from shapelab.common import moebius_inverse
  from shapelab.common import positive_angle_exponential as csqrt
  from shapelab.common.moebius_maps import specify_points as find_moebius
  from shapelab.common.moebius_maps import inverse_map as invert_moebius
end

if nargin<3
  interior = true(size(z));
  slit_interior = false(size(z));
  slit_exterior = slit_interior;
end
allz = interior | slit_interior | slit_exterior;

sgn = -sign(assumptions.winding_number);  % negative of winding number

z(interior) = csqrt(sgn*z(interior), 1/2);
vertices_int = real(csqrt(sgn*vertices_int, 1/2, 'cut_bias', sgn>0));
z(slit_interior) = real(csqrt(sgn*z(slit_interior), 1/2, 'cut_bias', sgn>0));
vertices_ext = real(csqrt(sgn*vertices_ext, 1/2, 'cut_bias', sgn<0));
z(slit_exterior) = real(csqrt(sgn*z(slit_exterior), 1/2, 'cut_bias', sgn<0));

% Find the terminal moebius map
% After the map we must have:
%    vertices_int(end-1) == -assumptions.tooth_length;
%    vertices_ext(end-1) == +assumptions.tooth_length;
%    vertices_int(end) == vertices_ext(end) ==0
% These three conditions uniquely specify a Moebius map
images = [0; -assumptions.tooth_length; assumptions.tooth_length];
assumptions.moebius_maps.terminal_map = invert_moebius(find_moebius([0; vertices_int(end-1); vertices_ext(end-1)], images));

% Now just apply the afore-calculated Moebius map
vertices_int = moebius_inverse(vertices_int, assumptions.moebius_maps.terminal_map);
vertices_ext = moebius_inverse(vertices_ext, assumptions.moebius_maps.terminal_map);
z(allz) = moebius_inverse(z(allz), assumptions.moebius_maps.terminal_map);
