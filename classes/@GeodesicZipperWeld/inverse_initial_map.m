function[z] = inverse_initial_map(self, z)
% inverse_initial_map -- The inverse initial map for zipper-type algorithms
%
% z = inverse_initial_map(self, z)
%
%     Performs the inverse of the initial map for zipper-type algorithms.

%z = moebius_inverse((-i*z).^2, mapdata.moebius_maps.initial_map);
z = self.moebius_maps.initial_map.inv((-i*z).^2);
