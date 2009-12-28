function[mapdata, z] = geodesic_driver(z, z_in, z_out, N_teeth, f_data, visdata);
% geodesic_driver -- The main driver for Geodesic algorithm computations
%
% mapdata = geodesic_driver(z, z_in, z_out, N_teeth, f_data, visdata)

persistent fa visualize
if isempty(fa)
  from shapelab.zipper.geodesic import base_map as fa
  from shapelab.zipper.drivers import visualize
end

z_in = [z_in; zeros([N_teeth 1])];
z_out = [z_out; zeros([N_teeth 1])];

a_array = zeros([N_teeth+1,1]);

interior_temp = zeros(size(z));
interior_temp(end) = 2;

boundary_temp = 2*ones([N_teeth+2 1]);
boundary_temp(end) = 1;

for q = 1:N_teeth
  a_array(q) = z(q);

  % Points still on the curve:
  f_data.point_id = interior_temp((q+1):end);
  z((q+1):end) = fa(z((q+1):end), a_array(q), f_data);

  % Points already unzipped:
  f_data.point_id = boundary_temp((N_teeth-q+2):end);
  [z_in(1:(q+1)),garbage] = fa(z_in(1:(q+1)),a_array(q),f_data);
  [garbage,z_out(1:(q+1))] = fa(z_out(1:(q+1)),a_array(q),f_data);

  if visdata.visualize
    visualize([0; z(q:(end-3))], z_in(1:(q+1)), z_out(1:(q+1)), visdata)
  end
end

z(1:N_teeth) = [];

[mapdata.a_array, mapdata.unzipped_in, mapdata.unzipped_out] = deal(...
         a_array,         z_in,                z_out);
