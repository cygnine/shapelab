function[z] = terminal_map(self, z)
% terminal_map -- The terminal map for zipper-type algorithms
%
% [z] = terminal_map(self, z)

%z = -sign(self.winding_number)*moebius(z, mapdata.moebius_maps.terminal_map).^2;
z = -sign(self.winding_number)*self.moebius_maps.terminal_map(z).^2;
