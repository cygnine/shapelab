function[z] = inverse_initial_map(z,mapdata)
% inverse_initial_map -- The inverse initial map for zipper-type algorithms
%
% [z] = inverse_initial_map(z, mapdata)
%
%     Performs the inverse of the initial map for zipper-type algorithms.

persistent moebius_inverse
if isempty(moebius_inverse)
  from shapelab.common import moebius_inverse
end

switch lower(mapdata.type)
case {'zipper', 'zipper_weld'}
  error('Not yet implemented');

  mapdata.moebius_maps.initial_map = [(z(2) - z(1))*[1, -z(3)]; ...
                              (z(2) - z(3))*[1, -z(1)]];

  z(1:3) = [Inf; -1; 0];
  z(4:end) = csqrt(moebius(z(4:end), mapdata.moebius_maps.initial_map));

case {'geodesic', 'slit', 'loewner'}
  z = moebius_inverse((-i*z).^2, mapdata.moebius_maps.initial_map);
otherwise
  error(['Unrecognized algorithm type: ' type]);
end
