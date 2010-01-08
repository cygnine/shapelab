function[w_n, z_n, mapdata] = terminal_map(w_n, z_n, mapdata)
% terminal_map -- The terminal map for zipper-type algorithms
%
% [w_n, z_n, mapdata] = terminal_map(w_n, z_n, mapdata)

persistent moebius
if isempty(moebius)
  from shapelab.common import moebius
end

switch lower(mapdata.type)
case {'geodesic', 'slit', 'zipper_weld'}
  mapdata.a_array(end) = 1/z_n(1);

  mapdata.terminal_map = [1                      0; ...
                          -mapdata.a_array(end)  1];

  w_n = -sign(mapdata.winding_number)*moebius(w_n, mapdata.terminal_map).^2;
  z_n = -sign(mapdata.winding_number)*moebius(z_n, mapdata.terminal_map).^2;
otherwise
  error('Unrecognized algorithm type');
end
