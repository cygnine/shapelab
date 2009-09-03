global handles;
clover = handles.shapelab.test_shapes.polar_clover;
gd = handles.shapelab.conformal_mapping;

N = 100;
opt.lobes = 4;
opt.lobe_depth_ratio = 0.95;
opt.lobe_max_radius = 1;
z = clover(N,opt);
zf = clover(10*N,opt);

opt.M = 50;
opt.node_locations = 'interior';
zint = clover(10*N,opt);

opt.node_locations = 'exterior';
zout = clover(10*N,opt);

mapdata = gd.compute_map_coordinates(z,'z_in', 0, 'w_in', 0, ...
      'type', 'zipper');

tin = unwrap(angle(mapdata.vertices_in));
tout = unwrap(angle(mapdata.vertices_out));

wint = zeros([N*M 1]);
for q = 1:M
  wint((q-1)*length(z)+1:q*length(z)) = exp(i*theta).*(q-1)/M;
end
wout = zeros([N*M 1]);
for q = 1:M
  wout((q-1)*length(z)+1:q*length(z)) = exp(i*theta).*(q+M)/M;
end
wint_image = gd.evaluate_inverse_map(wint,mapdata);
wout_image = gd.evaluate_inverse_map(wout,mapdata);
%
%wout = zeros([N*M 1]);
%for q = 1:M
%  wout((q-1)*length(z)+1:q*length(z)) = exp(i*theta).*(q+M)/M;
%end

unzipped_in_fine = exp(i*thetaf).';
unzipped_in_fine_image = gd.switch_zipper_side(unzipped_in_fine, mapdata, 'point_id',...
  ones(size(unzipped_in_fine)));
unzipped_in_fine_shape = gd.evaluate_inverse_map(unzipped_in_fine, mapdata,...
          'point_id', ones(size(unzipped_in_fine)));

unzipped_out_fine = exp(i*thetaf).';
unzipped_out_fine_image = gd.switch_zipper_side(unzipped_out_fine, mapdata, 'point_id',...
  2*ones(size(unzipped_out_fine)));
unzipped_out_fine_shape = gd.evaluate_inverse_map(unzipped_out_fine, mapdata,...
          'point_id', ones(size(unzipped_out_fine)));
