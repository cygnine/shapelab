function[w_n, z_n, mapdata] = moebius_alignment(w_n, z_n, mapdata)
% moebius_alignment -- Application of final Moebius map alignments for conformal maps
%
% [w_n, z_n, mapdata] = moebius_alignment(w_n, z_n, mapdata)

persistent moebius maps
if isempty(moebius)
  from shapelab.common import moebius

  maps.H_to_D = [-1, i; ...
                  1, i];
end

w_n = moebius(w_n, mapdata.moebius_maps.interior_terminal);
z_n = moebius(z_n, mapdata.moebius_maps.exterior_terminal);

% Finally, map to unit circle
w_n = moebius(w_n, maps.H_to_D);
z_n = moebius(z_n, maps.H_to_D);
