function[z] = backward_map(z,mapdata,varargin)
% backward_map -- Maps the disc to a shape
%
% [z] = backward_map(z,mapdata,{point_id=false(size(z))})
%
%     Uses the map mapdata to take points from the disc to the interior of the
%     shape. Assumes a welding map so it maps 'exterior' points as well. Whether
%     a point is 'exterior' or 'interior' is not needed as input, but the map
%     is different depending on this condition. 
%
%     point_id == false ---> points are 'interior' to the disc
%     point_id == true  ---> points are 'exterior' to the disc
%
%     Note that there is no ambiguity unless abs(z)==1. For those z such that
%     abs(z)~=1, the input point_id is ignored.
%
%     This function is the opposite of unzip_shape.
%
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%     conformal mapping", 2006.

persistent strict_inputs zip moebius 
persistent inverse_initial_map inverse_terminal_map assert_enough_points select_slider inverse_moebius_alignment

if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.common import moebius 
  from shapelab.zipper import inverse_initial_map inverse_terminal_map assert_enough_points
  from shapelab.zipper import select_slider inverse_moebius_alignment
end

zsize = size(z);
z = z(:);

opt = strict_inputs({'point_id'}, {false(size(z))}, [], varargin{:});

% Set point_id -- override user specifications if magnitude of z dictates
point_id = opt.point_id(:);
tol = 1e-8;
nflags = abs(z)<(1-tol);
point_id(nflags) = false;
pflags = abs(z)>(1+tol);
point_id(pflags) = true;

interior = zeros(size(z));
interior = (nflags | pflags);
slit_interior = not(interior) & not(point_id);
slit_exterior = not(interior) & point_id;

% We're doing all the maps backwards now:
[z(~point_id), z(point_id)] = inverse_moebius_alignment(z(~point_id), z(point_id), mapdata);

% Fix real-valued crap:
z(slit_interior) = real(z(slit_interior));
z(slit_exterior) = real(z(slit_exterior));

z = inverse_terminal_map(z, mapdata, interior, slit_interior, slit_exterior);

slide = select_slider(mapdata.type);

for q = mapdata.N_teeth:-1:1
  [z] = slide('zipup', q, z, mapdata, interior, slit_interior, slit_exterior);
end

z = reshape(inverse_initial_map(z, mapdata), zsize);
