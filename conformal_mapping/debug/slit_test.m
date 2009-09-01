global handles;
gd = handles.shapelab.conformal_mapping;

N = 100;
theta = linspace(-pi,pi,N+1); 
theta = theta(1:N);
z = (1+0.95*cos(theta*4)).*exp(i*theta);

thetaf = linspace(-pi,pi,10*N);
zf = (1+0.95*cos(thetaf*4)).*exp(i*thetaf);
M = 50;
zint = zeros([N*M,1]);
for q = 1:M
  zint((q-1)*length(z)+1:q*length(z)) = z.*((q)/(M+1));
end

zout = zeros([N*M,1]);
for q = 1:M
  zout((q-1)*length(z)+1:q*length(z)) = z.*((q+M+1)/M); 
end

zip_magnitude = 0.85;

mapdata = gd.compute_map_coordinates(z,'z_in', 0, 'w_in', 0, ...
      'zip_magnitude', zip_magnitude,'type', 'slit');

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
