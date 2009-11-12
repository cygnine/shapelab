from shapelab.test_shapes import polar_clover
from shapelab.common import polar_linspace as plinspace
imp shapelab.conformal_mapping.zipper as gd
imp shapelab.welding as weld

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

mapdata = gd.compute_map_coordinates(z,'z_in', 0.8, 'w_in', 0, ...
      'z_out', -1.1*i, 'w_out', Inf, 'type', 'zipper_weld');

tin = unwrap(angle(mapdata.vertices_in));
tout = unwrap(angle(mapdata.vertices_out));

wint = plinspace(opt.M, N, 'r0',0, 'r1', (opt.M-1)/opt.M);
wout = plinspace(opt.M, N, 'r0',1+1/opt.M, 'r1', 2);

wint_image = gd.evaluate_inverse_map(wint,mapdata);
wout_image = gd.evaluate_inverse_map(wout,mapdata);

unzipped_in_fine = plinspace(1,10*N,'r0',1,'r1',1);
unzipped_in_fine_image = weld.switch_zipper_side(unzipped_in_fine, mapdata, 'point_id',...
  ones(size(unzipped_in_fine)));
unzipped_in_fine_shape = gd.evaluate_inverse_map(unzipped_in_fine, mapdata,...
          'point_id', ones(size(unzipped_in_fine)));

unzipped_out_fine = plinspace(1,10*N,'r0',1,'r1',1);
unzipped_out_fine_image = weld.switch_zipper_side(unzipped_out_fine, mapdata, 'point_id',...
  2*ones(size(unzipped_out_fine)));
unzipped_out_fine_shape = gd.evaluate_inverse_map(unzipped_out_fine, mapdata,...
          'point_id', 2*ones(size(unzipped_out_fine)));
