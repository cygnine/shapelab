% This script picks out one shape from the mpeg7 shape database and attempts to
% approximate it. 

% From whichever directory you run this script, the mpeg7_contours.mat file
% should be present.

clear
close all;

global packages;
explot = packages.labtools.explot;  % Some common functions
ltex = packages.labtools.typelatex;

shapelab = packages.shapelab;  % shapelab stuff
zipper = shapelab.conformal_mapping.zipper;
welding = packages.shapelab.welding;
clover = shapelab.test_shapes.polar_clover;
fprint_norm = welding.normalize_fingerprint;

load mpeg7_contours;

shape = mpeg7_contour{3};
N = 150;

z = shape(round(linspace(1,length(shape),N+1)));
z(end) = [];

map_opt.type = 'geodesic';
map_opt.shape_0 = mean(shape);
fprintf('Computing geodesic map....\n');
geo_map = welding.make_zipper_map(z,map_opt);

map_opt.type = 'slit';
fprintf('Computing slit map....\n');
slit_map = welding.make_zipper_map(z,map_opt);

map_opt.type = 'zipper_weld';
fprintf('Computing zipper map....\n');
zip_map = welding.make_zipper_map(z,map_opt);

w_fine = shapelab.common.polar_linspace(1,100*N,'r0', 1, 'r1', 1);
point_id = ones(size(w_fine));  % interior of unit circle
z_geo_fine = zipper.evaluate_inverse_map(w_fine, geo_map,'point_id', point_id);
z_slit_fine = zipper.evaluate_inverse_map(w_fine, slit_map,'point_id', point_id);
z_zip_fine = zipper.evaluate_inverse_map(w_fine, zip_map,'point_id', point_id);
point_id = 2*ones(size(w_fine)); % exterior of unit circle
z_geo_fine = [z_geo_fine; zipper.evaluate_inverse_map(w_fine, geo_map,'point_id', point_id)];
z_slit_fine = [z_slit_fine; zipper.evaluate_inverse_map(w_fine, slit_map,'point_id', point_id)];
z_zip_fine = [z_zip_fine; zipper.evaluate_inverse_map(w_fine, zip_map,'point_id', point_id)];

figure(3);
subplot(2,2,1);
title('The geodesic approximation')
plot(z_geo_fine, 'r.'); hold on; plot(shape, 'b'); set(plot(z, 'k.'), 'markersize', 10);
subplot(2,2,2);
title('The slit approximation')
plot(z_slit_fine, 'r.'); hold on; plot(shape, 'b'); set(plot(z, 'k.'), 'markersize', 10);
subplot(2,2,3);
title('The zipper approximation')
plot(z_zip_fine, 'r.'); hold on; plot(shape, 'b'); set(plot(z, 'k.'), 'markersize', 10);
subplot(2,2,4);
plot(0,0, 'r.'); hold on; plot(0,0,'b'); hold on; plot(0,0,'k.');
axis([1,2,1,2]);
legend('Mapped inverse values', 'Original shape', 'Original data points');
