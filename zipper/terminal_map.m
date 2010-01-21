function[z] = terminal_map(z, mapdata)
% terminal_map -- The terminal map for zipper-type algorithms
%
% [z] = terminal_map(z, mapdata)

persistent moebius
if isempty(moebius)
  from shapelab.common import moebius
end

switch lower(mapdata.type)
case {'geodesic', 'slit', 'zipper_weld','loewner'}
  z = -sign(mapdata.winding_number)*moebius(z, mapdata.moebius_maps.terminal_map).^2;
otherwise
  error('Unrecognized algorithm type');
end
