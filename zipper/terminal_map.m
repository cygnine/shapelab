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
case 'loewner'
  % Same as above, but first rotate so that z_n(N) is at 0
  maps.H_to_D = [-1, i; ...
                  1, i];
  maps.D_to_H = [i, -i;...
                 -1, -1];

  rotation = -angle(moebius(z_n(end-3), maps.H_to_D));
  rotation = [exp(i*rotation), 0;...
              0,               1];

  map = maps.D_to_H*rotation*maps.H_to_D;

  w_n = moebius(w_n, map);
  z_n = moebius(z_n, map);

  % Now copy-paste from geodesic/slit/zipper_weld
  mapdata.a_array(end) = 1/z_n(1);

  mapdata.terminal_map = [1                      0; ...
                          -mapdata.a_array(end)  1];

  w_n = -sign(mapdata.winding_number)*moebius(w_n, mapdata.terminal_map).^2;
  z_n = -sign(mapdata.winding_number)*moebius(z_n, mapdata.terminal_map).^2;
otherwise
  error('Unrecognized algorithm type');
end
