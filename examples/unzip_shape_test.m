% This is a script to test unzipping algorithms

clear
close all

load('examples/mpeg7_contours.mat');

% This shape is relatively easy
samples = mpeg7_contour{5};

% As usual, slit messes up. Each of the following two screw slit up (for N=150).
% However, both geodesic and loewner power through.
%samples = mpeg7_contour{329};
%samples = mpeg7_contour{514};

N = 150;
z = samples(round(linspace(1,length(samples),N+1)));
z(end) = [];

from shapelab.zipper import unzip_shape

geo_map = unzip_shape(z, 'visualize', true);
slit_map = unzip_shape(z, 'type', 'slit', 'visualize', true);
loewner_map = unzip_shape(z, 'type', 'loewner', 'visualize', true);

figure;
plot(unwrap(angle(geo_map.vertices_in)), unwrap(angle(geo_map.vertices_out)), 'b.-'); hold on;
plot(unwrap(angle(slit_map.vertices_in)), unwrap(angle(slit_map.vertices_out)), 'r.-');
plot(unwrap(angle(loewner_map.vertices_in)), unwrap(angle(loewner_map.vertices_out)), 'k.-');

title('Welding map');
legend('Geodesic', 'Slit', 'Loewner');
