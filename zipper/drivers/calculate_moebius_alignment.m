function[w_n, z_n, mapdata] = calculate_moebius_alignment(w_n, z_n, mapdata)
% calculate_moebius_alignment -- Calculate and apply final Moebius map alignments for conformal maps
%
% [w_n, z_n, mapdata] = calculate_moebius_alignment(w_n, z_n, mapdata)
%
%     These final Moebius maps do four things:
%       1.) Connect user-specified shape-exterior points. (default: maps
%           infinity outside shape to infinity outside disc) This is a self-map
%           on H.
%           (mapdata.moebius_maps.exterior_terminal)
%       2.) Connect user-specified shape-interior points. (default: do nothing)
%           This is a self-map on H.
%           (mapdata.moebius_maps.interior_terminal)
%       3.) Rotates both the interior and exterior welding map points so that
%           the first sample point has angle 0. These are self-maps on D.
%           (Already built into maps 1 and 2, but separately accessible in
%            mapdata.moebius_maps.interior_rotation and 
%            mapdata.moebius_maps.exterior_rotation)
%       4.) Maps the result of all these points to the disc. (A map from H to
%           D; not included in mapdata)

persistent moebius maps wrap
if isempty(moebius)
  from shapelab.common import moebius
  from labtools import interval_wrap as wrap

  maps.H_to_D = [-1, i; ...
                  1, i];
  maps.invert_D = [0, 1; ...
                   1, 0];
  maps.D_to_H = [i, -i;...
                 -1, -1];
end

% Always specify the exterior map
map_to = 1/mapdata.w_out;

z_out_on_disc = moebius(z_n(end-2), maps.invert_D*maps.H_to_D);
maps.map_z_out_on_disc_to_0 = ...
     [abs(z_out_on_disc)/z_out_on_disc, -abs(z_out_on_disc); ...
      -conj(z_out_on_disc),             1];

if abs(map_to)<(10*eps)
  maps.map_0_on_disc_to_goal = eye(2);
else
  maps.map_0_on_disc_to_goal = ...
     [1,             abs(map_to); ...
      conj(map_to),  abs(map_to)/map_to];
end

% To see what it does, follow it from the end
bigmap = maps.map_0_on_disc_to_goal*maps.map_z_out_on_disc_to_0*...
         maps.invert_D*maps.H_to_D;

% moebius(z_n(1), bigmap) is the location of the first sample point on the
% boundary of the disc. We want to rotate it to have angle 0.
exterior_rotation = -angle(moebius(z_n(1), bigmap));

mapdata.moebius_maps.exterior_rotation = [exp(i*exterior_rotation), 0; ...
                                          0,                        1];

% Again, follow it from the end
mapdata.moebius_maps.exterior_terminal = ...
      maps.D_to_H*maps.invert_D*mapdata.moebius_maps.exterior_rotation*bigmap;

% The above maps the real line to the real line...therefore normalize it to get
% rid of imaginary (machine eps crap) stuff
% God what I wouldn't give for a *= operator
mapdata.moebius_maps.exterior_terminal = real(mapdata.moebius_maps.exterior_terminal/...
                                         max(mapdata.moebius_maps.exterior_terminal(:)));

z_in_on_disc = moebius(z_n(end-1), maps.H_to_D);

if isempty(mapdata.z_in) & isempty(mapdata.w_in);
  maps.map_z_in_on_disc_to_0 = eye(2);
else
  maps.map_z_in_on_disc_to_0 = ...
     [abs(z_in_on_disc)/z_in_on_disc, -abs(z_in_on_disc); ...
      -conj(z_in_on_disc),            1];
end

bigmap = maps.map_z_in_on_disc_to_0*maps.H_to_D;

interior_rotation = -angle(moebius(z_n(1), bigmap));
mapdata.moebius_maps.interior_rotation = ...
       [exp(i*interior_rotation), 0; ...
        0,                        1];

mapdata.moebius_maps.interior_terminal = ...
     maps.D_to_H*mapdata.moebius_maps.interior_rotation*bigmap;

mapdata.moebius_maps.interior_terminal = ...
  real(mapdata.moebius_maps.interior_terminal/...
   max(mapdata.moebius_maps.interior_terminal(:)));


w_n = moebius(w_n, mapdata.moebius_maps.interior_terminal);
z_n = moebius(z_n, mapdata.moebius_maps.exterior_terminal);

% Before mapping to unit circle, also save fingerprint values
w = w_n(1:end-3);
mapdata.fprint_int = wrap(atan2(2*w, 1-w.^2), [0, 2*pi]);
w = z_n(1:end-3);
mapdata.fprint_ext = wrap(atan2(2*w, 1-w.^2), [0, 2*pi]);

% Finally, map to unit circle
w_n = moebius(w_n, maps.H_to_D);
z_n = moebius(z_n, maps.H_to_D);

% Fix machine eps crap:
mapdata.fprint_int(1) = 0;
mapdata.fprint_ext(1) = 0;
w_n(1) = 1;
z_n(1) = 1;
