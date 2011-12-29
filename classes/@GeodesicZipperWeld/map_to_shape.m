function[z] = map_to_shape(self, z, varargin)
% map_to_shape -- Maps the disc back to a shape
%
% z = map_to_shape(self, z, [[ point_id=false(size(z)) ]] )
%
%     Uses the welding map instance, this takes points from the disc to the interior of the
%     shape. Assumes a welding map so it maps 'exterior' points as well. Whether
%     a point is 'exterior' or 'interior' is not needed as input, but the map
%     is different depending on this condition. 
%
%       point_id == false ---> points are 'interior' to the disc
%       point_id == true  ---> points are 'exterior' to the disc
%
%     Note that there is no ambiguity unless abs(z)==1. For those z such that
%     abs(z)~=1, the input point_id is ignored.
%
%     This function basically implements in the inverse of the constructor.
%
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%     conformal mapping", 2006.

zsize = size(z);
z = z(:);
if length(varargin) > 0
  point_id = varargin{1}(:);
else
  point_id = false(size(z));
end

% Set point_id -- override user specifications if magnitude of z dictates
tol = 1e-8;
nflags = abs(z)<(1-tol);
point_id(nflags) = false;
pflags = abs(z)>(1+tol);
point_id(pflags) = true;

interior = zeros(size(z));
interior = (nflags | pflags);
slit_interior = not(interior) & not(point_id);
slit_exterior = not(interior) & point_id;

% For those z's in slit_interior that lie between self.interior_vertices(end)
% and self.interior_vertices(1), we'll have to make a special exception in
% self.inverse_terminal_map because of multi-valued power functions. We'll do
% the same for slit_exterior as well.
z_angles = angle(z);
slit_interior_limbo = slit_interior & (  ...
                (z_angles > self.interior_vertices(end)) & ...
                (z_angles < self.interior_vertices(1))  ...
              );
slit_exterior_limbo = slit_exterior & (  ...
                (z_angles > self.exterior_vertices(end)) & ...
                (z_angles < self.exterior_vertices(1))  ...
              );

% We're doing all the maps backwards now:
[z(~point_id), z(point_id)] = self.inverse_moebius_alignment(z(~point_id), z(point_id));

% Fix real-valued crap:
z(slit_interior) = real(z(slit_interior));
z(slit_exterior) = real(z(slit_exterior));

z = self.inverse_terminal_map(z, interior, slit_interior, slit_exterior, slit_interior_limbo, slit_exterior_limbo);

for q = self.N_teeth:-1:1
  %[z] = slide('zipup', q, z, mapdata, interior, slit_interior, slit_exterior);
  [z] = self.slider('zipup', q, z, interior, slit_interior, slit_exterior);
end

z = reshape(self.inverse_initial_map(z), zsize);
