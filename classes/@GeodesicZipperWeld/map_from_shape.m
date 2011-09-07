function[z] = map_from_shape(self, z)
% map_from_shape -- Maps the shape to the disc
%
% w = map_from_shape(self, z)
%
%     Uses the welding map instance, this takes points from the shape to the
%     interior or exterior of the disc. This cannot reliably detect points on
%     the shape, so the output points w will be strictly inside or outside the
%     disc, to machine precision.

zshape = size(z);
z = z(:);

z = self.initial_map(z);

for q = 1:self.N_teeth
  z = self.slider('unzip', q, z, true(size(z)), false(size(z)), false(size(z)));
end

z = self.terminal_map(z);
%flags = abs(z)<=1;
z = self.moebius_alignment(z);
%[z(flags), z(~flags)] = self.moebius_alignment(z(flags), z(~flags));

z = reshape(z, zshape);
