% Just a script to test out interpolation of a fingerprint

clear
close all

load('examples/mpeg7_contours.mat');
from shapelab.test_shapes import downsample_points
from shapelab.zipper import unzip_shape interpolate_fingerprint rechart_fingerprint
from shapelab.zipper import unzip_shape interpolate_weld
from labtools import interval_wrap as wrap

samples = mpeg7_contour{34};
z_in = 350+175*i;
z_in = 100+150*i;
w_in = 0;

N = 150;
z = downsample_points(samples, N);

mapdata = unzip_shape(z, 'w_in', w_in, 'z_in', z_in);
interval = [0, 2*pi];

vertices_int = wrap(angle(mapdata.vertices_in), interval);
vertices_ext = wrap(angle(mapdata.vertices_out), interval);

% Remove after debugging:
theta_int = linspace(0, 2*pi, 1e4);
%theta_int = linspace(vertices_int(37), vertices_int(38), 10);
%theta_int(1) = [];
theta_ext = linspace(0, 2*pi, 1e4);
%theta_ext = linspace(vertices_ext(37), vertices_ext(38), 10);

[psi_out, psi_in] = interpolate_fingerprint(vertices_int, vertices_ext, 'theta_int', theta_int, 'theta_ext', theta_ext);


% Let's try recharting:
%int_chart = [0, 2*pi/3, 4*pi/3];
%ext_chart = [0, 2*pi/3, 4*pi/3];
%
%int_chart = [0, 1, 5];
%ext_chart = [0, 2.5, 5];
%
%[vint, vext] = rechart_fingerprint(vertices_int, vertices_ext, int_chart, ext_chart);
%
%[p_out, p_in] = interpolate_fingerprint(vint, vext, 'theta_int', theta_int, 'theta_ext', theta_ext);


plot(theta_int, psi_out, 'b.', psi_in, theta_ext, 'b.', vertices_int, vertices_ext, 'r.')
hold on;
%plot(theta_int, p_out, 'k.', p_in, theta_ext, 'k.', vint, vext, 'g.');

%plot(int_chart, ext_chart, 'mo');

%x = linspace(0, 2*pi, 100);
%plot(x,x);
