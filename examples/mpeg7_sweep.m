% This script performs a crude estimate of each mpeg7 shape and makes a crude
% calculation of the error.

% From whichever directory you run this script, the mpeg7_contours.mat file
% should be present.

clear
close all;

global handles;
explot = handles.common.explot;  % Some common functions
ltex = handles.common.typelatex;

shapelab = handles.shapelab;  % shapelab stuff
zipper = shapelab.conformal_mapping.zipper;
welding = handles.shapelab.welding;
clover = shapelab.test_shapes.polar_clover;
fprint_norm = welding.normalize_fingerprint;

load mpeg7_contours;

N_shapes = length(mpeg7_contour);
N = 100; % Number of data points for each shape

geo_errors = cell(size(mpeg7_contour));
%slit_errors = cell(size(mpeg7_contour));
%zip_errors = cell(size(mpeg7_contour));

for q = 1:N_shapes
  fprintf('Computing shape %d\n', q);
  shape = mpeg7_contour{q};
  z = shape(round(linspace(1,length(shape),N+1)));
  z(end) = [];

  map_opt.type = 'geodesic';
  map_opt.shape_0 = mean(shape);
  geo_map = welding.make_zipper_map(z,map_opt);

  %map_opt.type = 'slit';
  %slit_map = welding.make_zipper_map(z,map_opt);

  %map_opt.type = 'zipper_weld';
  %fprintf('Computing zipper map....\n');
  %zip_map = welding.make_zipper_map(z,map_opt);

  w_fine = shapelab.common.polar_linspace(1,100*N,'r0', 1, 'r1', 1);
  point_id = ones(size(w_fine));  % interior of unit circle
  z_geo_fine = zipper.evaluate_inverse_map(w_fine, geo_map,'point_id', point_id);
  %z_slit_fine = zipper.evaluate_inverse_map(w_fine, slit_map,'point_id', point_id);
  %z_zip_fine = zipper.evaluate_inverse_map(w_fine, zip_map,'point_id', point_id);
  point_id = 2*ones(size(w_fine)); % exterior of unit circle
  z_geo_fine = [z_geo_fine; zipper.evaluate_inverse_map(w_fine, geo_map,'point_id', point_id)];
  %z_slit_fine = [z_slit_fine; zipper.evaluate_inverse_map(w_fine, slit_map,'point_id', point_id)];
  %z_zip_fine = [z_zip_fine; zipper.evaluate_inverse_map(w_fine, zip_map,'point_id', point_id)];

  geo_error = 0*shape;
  %slit_error = 0*shape;
  %zip_error = 0*shape;
  
  for qq = 1:length(shape);
    point_radius = abs(shape(qq) - map_opt.shape_0);
    geo_error(qq) = min(abs(z_geo_fine - shape(qq)))/point_radius;
    %slit_error(qq) = min(abs(z_slit_fine - shape(qq)))/point_radius;
    %zip_error(qq) = max(abs(z_zip_fine - shape(qq)));
  end

  geo_errors{q} = geo_error;
  %slit_errors{q} = slit_error;
  % zip_errors{q} = zip_error;
end

trouble_shapes = [];
randomly_chosen_tolerance = 5e-2;
for q = 1:length(mpeg7_contour);
  if (max(geo_errors{q})>randomly_chosen_tolerance)
    trouble_shapes(end+1) = q;
  end
end
