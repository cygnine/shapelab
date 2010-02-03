function[z_interior, z_exterior] = inverse_moebius_alignment(z_interior, z_exterior,mapdata)
% inverse_moebius_alignment -- Inverts the terminal Moebius alignment
%
% [z_interior, z_exterior] = inverse_moebius_alignment(z_interior, z_exterior,mapdata)
%
%     See the function moebius_alignment for a list of what this function
%     inverts.

persistent moebius_inverse moebius maps
if isempty(moebius)
  from shapelab.common import moebius moebius_inverse
  maps.H_to_D = [-1, i; ...
                  1, i];
  maps.D_to_H = [i, -i;...
                 -1, -1];
end

do_interior = not(isempty(z_interior));
do_exterior = not(isempty(z_exterior));

if do_interior
  z_interior = moebius(z_interior, maps.D_to_H);
  z_interior = moebius_inverse(z_interior, mapdata.moebius_maps.interior_terminal);
end

if do_exterior
  z_exterior = moebius(z_exterior, maps.D_to_H);
  z_exterior = moebius_inverse(z_exterior, mapdata.moebius_maps.exterior_terminal);
end
