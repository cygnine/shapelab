% Just a script to test out if zipper interpolation and recharting commute

clear
close all

load('examples/mpeg7_contours.mat');
from shapelab.test_shapes import downsample_points
from shapelab.zipper import unzip_shape interpolate_fingerprint rechart_fingerprint
from shapelab.zipper import unzip_shape interpolate_weld apply_rechart
from labtools import interval_wrap as wrap

samples = mpeg7_contour{34};
z_in = 350+175*i;
%z_in = 100+150*i;
w_in = 0;

N = 150;
z = downsample_points(samples, N);

mapdata = unzip_shape(z, 'w_in', w_in, 'z_in', z_in);
interval = [0, 2*pi];

% recharting specification:
int_chart = [0, 2*pi/3, 4*pi/3];
ext_chart = [0, 2*pi/3, 4*pi/3];

vertices_int = wrap(angle(mapdata.vertices_in), interval);
vertices_ext = wrap(angle(mapdata.vertices_out), interval);

[vint, vext,M] = rechart_fingerprint(vertices_int, vertices_ext, ...
int_chart, ext_chart);

% Remove after debugging:
theta_int = linspace(0, 2*pi, 1e3);
theta_ext = linspace(0, 2*pi, 1e3);

[psi_out, psi_in] = interpolate_fingerprint(vertices_int, vertices_ext, 'theta_int', theta_int, 'theta_ext', theta_ext);
[theta_int, psi_out] = apply_rechart(theta_int, psi_out,M);
[psi_in, theta_ext] = apply_rechart(psi_in, theta_ext, M);
% psi_in and psi_out were first interpolated then recharted


[psi_out_rechart, psi_in_rechart] = interpolate_fingerprint(vint, vext, 'theta_int', theta_int, 'theta_ext', theta_ext);
% psi_in_rechart and psi_out_rechart were first recharted then interpolated

plot(theta_int, psi_out, 'b.', psi_in, theta_ext, 'r.');
