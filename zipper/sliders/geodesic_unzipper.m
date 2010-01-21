function[w_n, z_n, mapdata] = geodesic_unzipper(tooth, z_n, w_n, mapdata)
% geodesic_unzipper -- Unzips one tooth using the geodesic algorithm
%
% [w_n, z_n, mapdata] = geodesic_unzipper(tooth, z_n, w_n, mapdata)
%
%     w_n -- points on the interior
%     z_n -- points on the exterior

persistent unzip moebius moebius_inverse
if isempty(unzip)
  from shapelab.loewner.solutions import normal_linear_slit_unzip as unzip
  from shapelab.common import moebius moebius_inverse
end

if tooth==1 % Then initialize data array
  N = length(z_n) - 3;
  mapdata.a_array = zeros([N-1 1]);
end

a_id = tooth + 2;
a = z_n(a_id);
mapdata.a_array(tooth) = a;

b = abs(a)^2/real(a);  % The point on real axis that the next map sends to infinity
c = abs(a)^2/imag(a);  % The image of a after the map

% Use a moebius map m to send the circular arc 0 -- a to the perpendicular line 
% 0 -- mapdata.tooth_length/c
map = [mapdata.tooth_length/c 0; -1/b 1];
z_n = moebius(z_n, map);
w_n = moebius(w_n, map);

% Now unzip the perpendicular line segment
[temp, temp2, temp3] = unzip(i*mapdata.tooth_length, ...
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
