function[z] = initial_map(z,mapdata)
% initial_map -- The initial map for zipper-type algorithms
%
% [z] = initial_map(z, mapdata)
%
%     Performs the initial map for zipper-type algorithms.

persistent csqrt moebius
if isempty(csqrt)
  from shapelab.common import positive_angle_square_root as csqrt
  from shapelab.common import moebius 
end

switch type
case {'zipper', 'zipper_weld'}
  z = csqrt(moebius(z, mapdata.moebius_maps.initial_map));
case {'geodesic', 'slit', 'loewner'}
  z = i*sqrt(moebius(z, mapdata.moebius_maps.initial_map));
otherwise
  error(['Unrecognized algorithm type: ' type]);
end
