% A script to test the BA extension

clear
close all
this_is_a_run_to_save_figures = false;
save_location = '/home/akil/work/notes/shape/ba-extension/images';
shape_name = 'triangle3';

load('mpeg7_contours.mat');

% For zipper:
from shapelab.test_shapes import downsample_points
from shapelab.zipper import unzip_shape backward_map
from labtools import interval_wrap as wrap

% For the BA extension + beltrami
from shapelab.extensions import weld_primitive_setup as wsetup
from shapelab.extensions import ba_fingerprint_driver as ba
from shapelab.extensions import fingerprint_extend
from shapelab.extensions import beltrami_coefficient_driver as mu
from shapelab.common import H_to_D_cover as phi
from shapelab.common import D_to_H_cover as inv_phi

% For plot outputs:
from labtools import sized_figsave figsave
from labtools import typelatex as ltex

samples = mpeg7_contour{613};
%z_in = 120 + 75*i; % shape 68
%z_in = 340 + 250*i; % shape 158
%z_in = 140 + 200*i; % shape 1377 % for sample (hard)
%z_in = mean(samples); % shape 1108
%z_in = 150 + 125*i; % shape 1021 % for sample (medium)
%z_in = mean(samples); % shape 1023  % for moving around z_in
%z_in = 200 + 269.5*i; % shape 612  % for moving around z_in
  %z_in = 250 + 269.5*i; % shape 612  % location 2
  z_in = 300 + 269.5*i; % shape 612  % location 3
%z_in = mean(samples); % shape 952 % for sample (easy)
%z_in = mean(samples); % shape 1169 % for sample (super easy)
w_in = 0;

N = 250;
z_shape = downsample_points(samples, N);
mapdata = unzip_shape(z_shape, 'w_in', w_in, 'z_in', z_in);
interval = [0, 2*pi];
% Extracting fingerprint:
theta_int = wrap(angle(mapdata.vertices_in), interval);
theta_ext = wrap(angle(mapdata.vertices_out), interval);
% Machine eps crap (not sure if this is necessary)
theta_int(1) = 0;
theta_ext(1) = 0;

% This just does precomputations to save time -- there are codes that do this
% already if you just input theta_int, theta_ext, but running these codes does
% certain (kinda expensive) computations repeatedly, which is wasteful.
stuff = wsetup(theta_int, theta_ext);

% Create mesh of points to map:
N_r = 20;
N_arg = 200;
%r = logspace(-3, -1e-5, N_r);
r = linspace(1e-3, 1-1e-3, N_r);
arg = linspace(0, 2*pi, N_arg+1);
w = (1-r)'*exp(i*arg);

circle = exp(i*arg);

% Plots: 
% 1.) conformal maps
% 2.) fingerprint
% 3.) BA extension (disc)
% 4.) Beltrami coefficient

% Conformal map images:
conformal_interior = backward_map(w, mapdata);
conformal_exterior = backward_map(1./w, mapdata, 'point_id', true(size(w)));

% Plotting:
conformal_plot = figure;
plot(z_shape, 'r.-'); hold on;
plot(conformal_interior, 'k'); plot(conformal_interior.', 'k'); 
temp = axis;
axis equal;
plot(conformal_exterior, 'k'); plot(conformal_exterior.', 'k');
temp = [2*(temp(1:2)-mean(temp(1:2))) + mean(temp(1:2)), ...
        2*(temp(3:4)-mean(temp(3:4))) + mean(temp(3:4))]; axis(temp); axis off;

% A new mesh for a new plot:
N_r = 40;
N_arg = 300;
% Points on the disc:
r = logspace(-8, -1e-5, N_r);
%r = linspace(1e-8, 1-1e-3, N_r);
arg = linspace(0, 2*pi, N_arg+1);
w = (1-r)'*exp(i*arg);
z = inv_phi(w);

% Interior extension
H_int = phi(ba(stuff, z));

% New mesh
N_r = 20;
N_arg = 1200;
r = linspace(1e-3, 1-1e-3, N_r);
arg = linspace(0, 2*pi, N_arg+1);
w = (1-r)'*exp(i*arg);
z = inv_phi(w);

% Exterior extension
H_ext = phi(ba(stuff, conj(z)));

% Plotting
extension_plot = figure; plot(H_int, 'k'); hold on; plot(H_int.', 'k'); 
plot(circle, 'r', 'linewidth', 3)
axis equal; 
temp = axis;
plot(H_ext, 'k'); plot(H_ext.', 'k');
axis([-2 2 -2 2]);
axis off;

% Plotting fingerprint
weld_plot = figure;
x = linspace(-pi, 3*pi, 1e3);
plot(x, fingerprint_extend(stuff,x), 'k', 'linewidth', 3);
temp = axis;
axis([-pi 3*pi temp(3:4)]);
ltex(xlabel('$x$', 'fontsize', 20, 'fontweight', 'b'));
ltex(ylabel('$h(x)$', 'fontsize', 20, 'fontweight', 'b'));

% A pretty fine mesh is needed to get rid of artifacts for badly-behaved
% beltrami coefficients
N_r = 700;
N_arg = 700;
r = linspace(1e-3, 1-1e-3, N_r);
arg = linspace(0, 2*pi, N_arg+1);
w = (1-r)'*exp(i*arg);

% Compute beltrami coefficient
m = mu(stuff, w);

% Plotting:
mu_plot = figure;
surfc(real(w), imag(w), abs(m), 'EdgeColor', 'none'); view([0 0 1]); axis equal; axis off

% Saving data:
if this_is_a_run_to_save_figures 
  plot_names = {'conformal', 'extension', 'weld', 'mu'};
  for q = 1:4
    plot_names{q} = strcat(shape_name, '_', plot_names{q});
  end
  sized_figsave([conformal_plot, extension_plot, mu_plot], ...
                plot_names([1,2,4]), save_location);
  figsave([weld_plot], plot_names(3),save_location);
end
