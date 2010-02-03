function[z, mapdata] = calculate_initial_map(z,type,mapdata)
% calculate_initial_map -- The initial map for zipper-type algorithms
%
% [z, mapdata] = calculate_initial_map(z, type, mapdata)
%
%     Determines and computes the initial map for zipper-type algorithms.

persistent csqrt moebius
if isempty(csqrt)
  from shapelab.common import positive_angle_square_root as csqrt
  from shapelab.common import moebius 
end

switch type
case {'zipper', 'zipper_weld'}
  mapdata.moebius_maps.initial_map = [(z(2) - z(1))*[1, -z(3)]; ...
                              (z(2) - z(3))*[1, -z(1)]];

  z(1:3) = [Inf; -1; 0];
  z(4:end) = csqrt(moebius(z(4:end), mapdata.moebius_maps.initial_map));

case {'geodesic', 'slit', 'loewner'}
  mapdata.moebius_maps.initial_map = [1 -z(2); ...
                              1 -z(1)];
  z(1:2) = [Inf; 0];
  z(3:end) = i*sqrt(moebius(z(3:end), mapdata.moebius_maps.initial_map));

otherwise
  error(['Unrecognized algorithm type: ' type]);
end
