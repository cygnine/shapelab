function[z] = geodesic_slider(direction, tooth, z, mapdata, interior, slit_interior, slit_exterior)
% geodesic_slider -- Unzips one tooth using the geodesic algorithm
%
% [z] = geodesic_slider(direction, tooth, z, mapdata, ...
%          interior, slit_interior, slit_exterior)
%
%     The three final inputs are mutually exclusive boolean flags that
%     determine how each point z is treated:
%
%       interior==true ----> z is somewhere in the interior, no special treatment
%                            is necessary.
%       slit_interior==true ----> z lies on the slit, on the interior part
%       slit_exterior==true ----> z lies on the slit, on the exterior part
%
%       direction -- either 'unzip' or 'zipup'
%
%     If the three final inputs are omitted, the default is to treat all points
%     as interior points.

persistent unzip zipup moebius moebius_inverse
if isempty(unzip)
  from shapelab.loewner.solutions import normal_linear_slit_unzip_wrapper as unzip
  from shapelab.loewner.solutions import normal_linear_slit_zip as zipup
  from shapelab.common import moebius moebius_inverse
end

switch direction
case 'unzip'
  if nargin<5
    interior = true(size(z));
    slit_exterior = ~interior;
    slit_interior = slit_exterior;
  end

  a = mapdata.a_array(tooth);
  b = abs(a)^2/real(a);  
  c = abs(a)^2/imag(a); 
  map = [mapdata.tooth_length/c 0; -1/b 1];

  z = moebius(z, map);

  [z(interior), z(slit_interior), z(slit_exterior)] = ...
    unzip(i*mapdata.tooth_length, z(interior), z(slit_interior), z(slit_exterior));

case 'zipup'
  a = mapdata.a_array(tooth);
  b = abs(a)^2/real(a);
  c = abs(a)^2/imag(a);  
  map = [mapdata.tooth_length/c 0; -1/b 1];

  z = zipup(i*mapdata.tooth_length,z); 

  z = moebius_inverse(z, map);
end
