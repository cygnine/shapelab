% This is a test function for computing a shape welding map using a Loewner
% evolution.

clear
close all
load mpeg7_contours
from shapelab.zipper import construct_map

samples = mpeg7_contour{73};
N = 400;

z = samples(round(linspace(1,length(samples),N+1)));
z(end) = [];

geodesic_map = construct_map(z, 'type', 'geodesic', 'visualize', true);
loewner_map = construct_map(z, 'type', 'loewner', 'visualize', true);

% Test to see if the welding map is accurate
temp = diff(unwrap(angle(loewner_map.vertices_in)));
if any(temp<0)
  fprintf('Some vertices_in are decreasing\n');
end
temp = diff(unwrap(angle(loewner_map.vertices_out)));
if any(temp<0)
  fprintf('Some vertices_out are decreasing\n');
end

figure; plot(unwrap(angle(geodesic_map.vertices_in)), unwrap(angle(geodesic_map.vertices_out)), 'r.-', ...
                unwrap(angle(loewner_map.vertices_in)), unwrap(angle(loewner_map.vertices_out)), 'b.-')
