% This script file tests various capabilities of the zipper transformation:
%
% - computation of a welding map
% - inverting the map (disc ---> shape)
% - interpolating a welding map given the mapdata
% - interpolating a welding map without the mapdata

clear 
close all

load('examples/mpeg7_contours.mat');
from shapelab.test_shapes import downsample_points
from shapelab.zipper import unzip_shape backward_map
from shapelab.common import polar_linspace

samples = mpeg7_contour{54};
% For this shape:
z_in = 100+150*i;
z_in = 350+100*i;
z_in = 360+20*i;
w_in = 0;

N = 150;
z = downsample_points(samples, N);
z = circshift(z, 102);

mapdata = unzip_shape(z, 'w_in', w_in, 'z_in', z_in);

w = polar_linspace(40, 40, 'r0', 0.01, 'r1', 1);
wz = backward_map(w, mapdata);

w_circ = exp(i*linspace(0, 2*pi, 1e4).');
wz_circ = backward_map(w_circ, mapdata, 'point_id', false(size(w_circ)));

wv = polar_linspace(40, 80, 'r0', 1, 'r1', 3);
wvz = backward_map(wv, mapdata, 'point_id', true(size(wv)));

plot(wz, 'r-'); hold on; plot(wvz, 'r-'); plot(wz_circ, 'g.-'); plot(z, 'b.')
plot(z(1), 'kx');
