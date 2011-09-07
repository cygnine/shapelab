function[z_int, z_ext] = inverse_moebius_alignment(self, z_int, z_ext)
% inverse_moebius_alignment -- Inverts the terminal Moebius maps
%
% [z_int, z_ext] = inverse_moebius_alignment(self, z_int, z_ext)

%do_interior = not(isempty(z_int));
%do_exterior = not(isempty(z_ext));

%if do_interior
%  z_interior = moebius(z_interior, maps.D_to_H);
%  z_interior = moebius_inverse(z_interior, mapdata.moebius_maps.interior_terminal);
%end
%if do_exterior
%  z_exterior = moebius(z_exterior, maps.D_to_H);
%  z_exterior = moebius_inverse(z_exterior, mapdata.moebius_maps.exterior_terminal);
%end

z_int = self.moebius_maps.interior_rotation.inv(z_int);
z_int = self.moebius_maps.D_to_H(z_int);
z_int = self.moebius_maps.interior_terminal.inv(z_int);

z_ext = self.moebius_maps.exterior_rotation.inv(z_ext);
z_ext = self.moebius_maps.D_to_H(z_ext);
z_ext = self.moebius_maps.exterior_terminal.inv(z_ext);
