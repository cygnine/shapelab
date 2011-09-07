function[w, z] = moebius_alignment(self, w, z)
% moebius_alignment -- Application of final Moebius map alignments for conformal maps
%
% [w, z] = moebius_alignment(self, w, z)
%
%     w are interior points
%     z are exterior points

if nargin < 3
  w2 = self.moebius_maps.interior_terminal(w);
  w2 = self.moebius_maps.H_to_D(w2);
  flags = abs(w2)>=1;
  w2 = self.moebius_maps.interior_rotation(w2);

  % Cheap and dumb:
  w2(flags) = self.moebius_maps.exterior_terminal(w(flags));
  w2(flags) = self.moebius_maps.H_to_D(w2(flags));
  w2(flags) = self.moebius_maps.exterior_rotation(w2(flags));

  w = w2;
else
  w = self.moebius_maps.interior_terminal(w);
  w = self.moebius_maps.H_to_D(w);
  w = self.moebius_maps.interior_rotation(w);

  z = self.moebius_maps.exterior_terminal(z);
  z = self.moebius_maps.H_to_D(z);
  z = self.moebius_maps.exterior_rotation(z);
end
