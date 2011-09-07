function[mapdata] = unzip_shape(z_n,varargin)
% unzip_shape -- `Unzips' a shape via a zipper-type algorithm
%
% [mapdata] = unzip_shape(z_n,{z_in=false,w_in=false,
%                 z_out=false, w_out=false, winding_number=1,
%                 zip_magnitude=0.85, type='geodesic', M=30})
%
%     The main workhorse routine for computing a welding map of a shape.
%
%     Output mapdata contains:
%            a_array
%            moebius_maps
%            z_in
%            w_in
%            z_out
%            w_out
%            vertices_in
%            vertices_out
%            type
%
%     The output mapdata contains all the necessary information for travelling
%     back and forth on the conformal map defined by the shape's vertices z_n
%     and their corresponding inner and outer values on the unit circle,
%     vertices_in and vertices_out.
%
%     This is the main work file for zipper-type algorithms from [1]. The 'base'
%     map may be interchanged to form the 'geodesic', 'slit', and 'zipper'
%     algorithms.
%    
%     The map takes z_in (inside the original shape) to w_in (inside the unit
%     circle). If z_in is specified without w_in, then w_in is assumed to be 0.
%     If neither z_in nor w_in is specified, no interior Moebius transform is
%     performed. 
%
%     The map takes z_out (outside the original shape) to w_out (outside the
%     unit circle.) If z_out is specified without w_out, then w_out is assumed
%     to be infinity. If neither z_out nor w_out is specified, the map takes
%     infinity to infinity.
%
%     The winding_number determines the orientation of the input data with
%     respect to the shape.
%
%     The option zip_magnitude refers to how each 'base' map is normalized.
%     Empirically, values O(1) produce stable results.
%
%     The option M is the number of sub-steps that a predictive Loewner solver
%     takes between each shape vertex.
%
%     The `type' of mapping can be: 'geodesic', 'slit', 'zipper', or 'loewner'.
%
%     A note on normalization: this function normalizes the interior and
%     exterior moebius maps so that the first vertex of the shape lies at the
%     mapped point z=1 (both interior and exterior).
%
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%          conformal mapping", 2006.

persistent input_schema zip moebius visualize
persistent initial_map terminal_map assert_enough_points select_unzipper moebius_alignment

if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import moebius 
  from shapelab.zipper import assert_enough_points
  from shapelab.zipper import select_unzipper
  from shapelab.zipper.drivers import visualize
  from shapelab.zipper.drivers import calculate_initial_map as initial_map
  from shapelab.zipper.drivers import calculate_terminal_map as terminal_map
  from shapelab.zipper.drivers import calculate_moebius_alignment as moebius_alignment
end

inputs = {'z_in', 'w_in', 'z_out', 'w_out', 'winding_number',...
          'zip_magnitude','type','visualize', 'M'};
defaults = {[], [], Inf, Inf, 1, 0.85, 'geodesic',false, 30};
opt = input_schema(inputs,defaults,[],varargin{:});

mapdata.w_out = opt.w_out;
mapdata.w_in = opt.w_in;
mapdata.z_in = opt.z_in;
mapdata.z_out = opt.z_out;
mapdata.winding_number = opt.winding_number;
mapdata.type = opt.type;
mapdata.tooth_length = opt.zip_magnitude;
mapdata.M = opt.M;

moebius_maps.H_to_D = [-1, i; ...
                        1, i];  % Just a map from half plane to unit circle.

% Some initial computations, making sure we have enough points, and append some
% more points
if isempty(opt.z_in)
  [N, N_teeth, z_n] = assert_enough_points(z_n, opt.type, opt.z_out, mean(z_n), z_n(1));
else
  [N, N_teeth, z_n] = assert_enough_points(z_n, opt.type, opt.z_out, opt.z_in, z_n(1));
end

% Now length(z_n) == N+3
[z_n, mapdata] = initial_map(z_n, opt.type, mapdata);

w_n = z_n;
% z_n are the 'exterior' points, and w_n are the interior points

visdata = initialize_visualization();
visdata.visualize = opt.visualize;

%% Initialization for looping over teeth
unzip = select_unzipper(opt.type);

% Loop over teeth:
for q = 1:N_teeth
  [w_n, z_n, mapdata] = unzip(q, z_n, w_n, mapdata);
  visualize(z_n(q+2:end-3), w_n(1:q+2), z_n(1:q+2), visdata);
end

[w_n, z_n, mapdata] = terminal_map(w_n, z_n, mapdata);

[w_n, z_n, mapdata] = moebius_alignment(w_n, z_n, mapdata);

% Done
mapdata.vertices_in = w_n(1:end-3);
mapdata.vertices_out = z_n(1:end-3);

  function[visdata]=initialize_visualization()
    if opt.visualize
      figure(); 
      theta = linspace(0,2*pi,200);
      plot(exp(i*theta), 'b--');
      axis square;
      axis off;
      hold on;
      shape_plot = plot(moebius(z_n, moebius_maps.H_to_D), 'k.-');
      temp = moebius([Inf; 0],moebius_maps.H_to_D);
      zin_plot = plot(real(temp), imag(temp), 'r.');
      temp = moebius([Inf; 0], moebius_maps.H_to_D);
      zout_plot = plot(real(temp), imag(temp), 'b.');
    else
      shape_plot = 0;
      zin_plot = 0;
      zout_plot = 0;
    end

    visdata.shape_plot = shape_plot;
    visdata.zin_plot = zin_plot;
    visdata.zout_plot = zout_plot;
  end

end
