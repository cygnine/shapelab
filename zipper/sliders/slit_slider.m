function[w_n, z_n, mapdata] = slit_slider(direction, tooth, z_n, w_n, mapdata)
% slit_slider -- Unzips one tooth using the slit algorithm
%
% [w_n, z_n, mapdata] = slit_slider(direction, tooth, z_n, w_n, mapdata)

persistent unzip zipup 
if isempty(unzip)
  from shapelab.loewner.solutions import oblique_linear_slit_unzip as unzip
  from shapelab.loewner.solutions import oblique_linear_slit_zip as zipup
end

switch direction
case 'down'
  if tooth==1 % Then initialize data array
    N = length(z_n) - 3;
    mapdata.a_array = zeros([N-1 1]);
    mapdata.a_angle = zeros([N-1 1]);
  end

  a_id = tooth + 2;
  a = z_n(a_id);
  mapdata.a_array(tooth) = a;
  mapdata.a_angle(tooth) = angle(a);

  [temp, temp2, temp3] = unzip(a, angle(a), ...
                               [z_n(1:a_id-2); w_n(1:a_id-2); z_n(a_id+1:end)], ...
                               [z_n(a_id-1:a_id); w_n(a_id-1:a_id)]);

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
  error('Not yet implemented');
end
