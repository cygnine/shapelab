function[z] = slit_slider(direction, tooth, z, mapdata, interior, slit_interior, slit_exterior)
% slit_slider -- The discretized elementary conformal map using the slit algorithm
%
% [z] = slit_slider(direction, tooth, z, mapdata, ...
%         interior, slit_interior, slit_exterior)

persistent unzip zipup 
if isempty(unzip)
  from shapelab.loewner.solutions import oblique_linear_slit_unzip as unzip
  from shapelab.loewner.solutions import oblique_linear_slit_zip as zipup
end

switch direction
case 'unzip'
  if nargin<5
    interior = true(size(z));
    slit_interior = ~interior;
    slit_exterior = slit_interior;
  end

  a = mapdata.a_array(tooth);

  [z(interior), z(slit_interior), garbage] = ...
    unzip(a, angle(a), z(interior), z(slit_interior));

  [garbage, garbage, z(slit_exterior)] = ...
    unzip(a, angle(a), [], z(slit_exterior));

 % interior pts , 'left' points , 'right' points

  z_n(1:a_id-2) = temp(1:a_id-2);
  w_n(1:a_id-2) = temp(a_id-1:2*a_id-4);

  z_n(a_id+1:end) = temp(2*a_id-3:end);
  w_n(a_id+1:end) = temp(2*a_id-3:end);

  z_n(a_id-1:a_id) = temp3(1:2);
  w_n(a_id-1:a_id) = temp2(3:4);

  z_n(a_id) = 0;
  w_n(a_id) = 0;

case 'up'
  a = mapdata.a_array(tooth);
  temp = zipup(a, angle(a), [z_n; w_n]);
  z_n = temp(1:end/2);
  w_n = temp((end/2+1):end);
end
