% This script performs a crude estimate of each mpeg7 shape and makes a crude
% calculation of the error.

% From whichever directory you run this script, the mpeg7_contours.mat file
% should be present.

clear
close all;

global handles;
explot = handles.common.explot;  % Some common functions
ltex = handles.common.typelatex;
figsave = handles.common.figsave;

shapelab = handles.shapelab;  % shapelab stuff
zipper = shapelab.conformal_mapping.zipper;
welding = handles.shapelab.welding;
clover = shapelab.test_shapes.polar_clover;
fprint_norm = welding.normalize_fingerprint;

load mpeg7_contours;
%mpeg7_contour = mpeg7_contour(101:500);

N_shapes = length(mpeg7_contour);
N = 100; % Number of data points for each shape

geo_errors = cell(size(mpeg7_contour));
slit_errors = cell(size(mpeg7_contour));
zip_errors = cell(size(mpeg7_contour));

slit_map_errors = cell(0);
zip_map_errors = cell(0);

slit_calc = false;
zip_calc = false;

for q = 1:N_shapes
  fprintf('Computing shape %d\n', q);
  shape = mpeg7_contour{q};
  z = shape(round(linspace(1,length(shape),N+1)));
  z(end) = [];
  % Remove common nodes
  [z,zipable] = shapelab.common.shape_preprocessing(z);
  if not(zipable) % Then add a few points and try again
    z = shape(round(linspace(1,length(shape),N+11)));
    z(end) = [];
    [z,zipable] = shapelab.common.shape_preprocessing(z);
  end

  if zipable
    % If odd # of nodes, tack shape(2) into place -- in practice, this clause is
    % only needed for q = 1300
    if mod(length(z),2)==1
      z = [z(1); shape(2); z(2:end)];
    end

    map_opt.type = 'geodesic';
    map_opt.shape_0 = mean(shape);
    map_opt.visualize = false;
    fprintf('...geodesic...');
    geo_map = welding.make_zipper_map(z,map_opt);

    map_opt.type = 'slit';
    try
      fprintf('...slit...');
      slit_map = welding.make_zipper_map(z,map_opt);
      slit_calc = true;
    catch
      fprintf('...fail...');
      slit_map_errors{end+1} = q;
    end

    map_opt.type = 'zipper_weld';
    try
      fprintf('...zipper...');
      zip_map = welding.make_zipper_map(z,map_opt);
      zip_calc = true;
    catch
      fprintf('...FAIL...');
      zip_map_errors{end+1} = q;
    end
    fprintf('\n');

    fhandle = figure();
    subplot(2,2,1);
    w_fine = shapelab.common.polar_linspace(1,100*N,'r0', 1, 'r1', 1);
    point_id = ones(size(w_fine));  % interior of unit circle
    %z_geo_fine = zipper.evaluate_inverse_map(w_fine, geo_map,'point_id', point_id);
    z_geo_fine = [];
    point_id = 2*ones(size(w_fine)); % exterior of unit circle
    z_geo_fine = [z_geo_fine; zipper.evaluate_inverse_map(w_fine, geo_map,'point_id', point_id)];
    plot(z,'b.'); hold on; plot(z_geo_fine, 'r'); title('Geodesic');
    if slit_calc
      point_id = ones(size(w_fine));  % interior of unit circle
      %z_slit_fine = zipper.evaluate_inverse_map(w_fine, slit_map,'point_id', point_id);
      z_slit_fine = [];
      point_id = 2*ones(size(w_fine)); % exterior of unit circle
      z_slit_fine = [z_slit_fine; zipper.evaluate_inverse_map(w_fine, slit_map,'point_id', point_id)];
      subplot(2,2,2);
      plot(z,'b.'); hold on; plot(z_slit_fine, 'r'); title('Slit');
    end
    if zip_calc
      point_id = ones(size(w_fine));  % interior of unit circle
      %z_zip_fine = zipper.evaluate_inverse_map(w_fine, zip_map,'point_id', point_id);
      z_zip_fine = [];
      point_id = 2*ones(size(w_fine)); % exterior of unit circle
      z_zip_fine = [z_zip_fine; zipper.evaluate_inverse_map(w_fine, zip_map,'point_id', point_id)];
      subplot(2,2,3);
      plot(z,'b.'); hold on; plot(z_zip_fine, 'r'); title('Zipper');
    end

    figsave(fhandle,{num2str(q)}, 'mpeg7_data');

    geo_error = 0*shape;
    slit_error = 0*shape;
    zip_error = 0*shape;
    
    for qq = 1:length(shape);
      point_radius = abs(shape(qq) - map_opt.shape_0);
      geo_error(qq) = min(abs(z_geo_fine - shape(qq)))/point_radius;
      if slit_calc
        slit_error(qq) = min(abs(z_slit_fine - shape(qq)))/point_radius;
      end
      if zip_calc
        zip_error(qq) = max(abs(z_zip_fine - shape(qq)));
      end
    end

    geo_errors{q} = geo_error;
    if slit_calc
      slit_errors{q} = slit_error;
    end
    if zip_calc
      zip_errors{q} = zip_error;
    end
    
    slit_calc = false;
    zip_calc = false;
    close all
  end
end

%trouble_shapes = [];
%randomly_chosen_tolerance = 5e-2;
%for q = 1:length(mpeg7_contour);
%  if (max(geo_errors{q})>randomly_chosen_tolerance)
%    trouble_shapes(end+1) = q;
%  end
%end
