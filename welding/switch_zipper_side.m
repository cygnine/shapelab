function[z] = switch_zipper_side(z,mapdata,varargin)
% [z]= switch_zipper_side(z,mapdata,{point_id=ones(size(w))})
%
%     'Switches' the side of the zipper that the points z lie on. I.e., in
%     geometric language, this function zips ups the boundary points z according
%     to the map mapdata. It then unzips them again, but switches the side of
%     the zipper that the points lie on. In practice, this maps the unit circle
%     to the unit circle and is the fingerprint of the map. The optional input
%     point_id specifies the starting side of the points z.
%
%     z(point_id==1) -----> (default) points on the interior boundary of the
%     unit disc that you want to map onto the exterior.
%
%     z(point_id==2) -----> points on the exterior boundary of the unit disc
%     that you want to map onto the interior.

global handles;
opt = handles.common.InputSchema({'point_id'}, {ones(size(z))}, [], varargin{:});
shapelab = handles.shapelab;

switch lower(mapdata.type)
case 'geodesic'
  ifa = shapelab.conformal_mapping.geodesic.inverse_base_conformal_map;
  fa = shapelab.conformal_mapping.geodesic.base_conformal_map;
  zipper = false;
case 'slit'
  ifa = shapelab.conformal_mapping.slit.inverse_base_conformal_map;
  fa = shapelab.conformal_mapping.slit.base_conformal_map;
  zipper = false;
case 'zipper'
  error('This isn''t a welding-capable map');
case 'zipper_weld'
  ifa_geo = shapelab.conformal_mapping.geodesic.inverse_base_conformal_map;
  fa_geo = shapelab.conformal_mapping.geodesic.base_conformal_map;
  ifa = shapelab.conformal_mapping.zipper.inverse_base_conformal_map;
  fa = shapelab.conformal_mapping.zipper.base_conformal_map;
  zipper = true;
otherwise
  error(['Unrecognized map type ' mapdata.type]);
end

moebius = handles.shapelab.common.moebius;
moebius_inv = handles.shapelab.common.moebius_inverse;
csqrt = handles.shapelab.common.positive_angle_square_root;

N = length(mapdata.a_array)+1;
% We must deal with two cases:
zint_id = opt.point_id==1;
zext_id = opt.point_id==2;
z_tooth_indicator = zeros(size(z));
% A faster way to do the following would be to sort the input, store only the
% starting values of each bin, and then apply the inverse of the permutation
% operator at the end to return the sorted values to their unsorted initial
% state.
z_tooth_indicator_bins = cell([mapdata.N_teeth+1,1]);

if not(all(abs(abs(z)-1)<1e-10))
  error('It doesn''t look like you gave me points on the unit circle....');
end
% map unit circle to R:
m = [-1, i;...
      1, i];
z = real(moebius_inv(z,m));

% apply (inverse) normalizing moebius maps
z(zint_id) = moebius_inv(z(zint_id), mapdata.m_in);
z(zext_id) = moebius_inv(z(zext_id), mapdata.m_out);

% This slit/geodesic terminal map is the same for all cases
z = z*-sign(mapdata.winding_number);
z_copy = z;
z(z_copy<0) = i*sqrt(abs(z_copy(z_copy<0)));
temp = (z_copy>=0);
z(temp) = sign(opt.point_id(temp)-1.5).*sqrt(z_copy(temp));  % How's that for hackery?

% Invert the terminal map (moebius)
m = [1, 0;...
     -mapdata.a_array(end), 1];
z = moebius_inv(z,m);

ifa_opt.point_id = ones(size(z));
ifa_opt.point_id(abs(imag(z))>1e-10) = 0;
ifa_opt.cut_magnitude = mapdata.zip_magnitude;
opt_input = ifa_opt;

if zipper % Then do 1 inverse geodesic map
  unzipped = ifa_opt.point_id==1;
  opt_input.point_id = ifa_opt.point_id(unzipped);
  z_temp = ifa_geo(z(unzipped), mapdata.c_array(end), opt_input);

  % These points changed from being unzipped to being zipped_up
  newly_zipped = abs(imag(z_temp))>0;
  now_zipped = unzipped;
  temp = unzipped(unzipped);
  temp(not(newly_zipped)) = 0;
  now_zipped(unzipped) = temp;

  ifa_opt.point_id(now_zipped) = 0;
  z_tooth_indicator(now_zipped) = mapdata.N_teeth + 1;
  z_tooth_indicator_bins{mapdata.N_teeth+1} = find(now_zipped);
  z(unzipped) = z_temp;
end

% Zip up all the zipper teeth
for q = mapdata.N_teeth:-1:1
  % Only invert the map for points that are unzipped
  unzipped = ifa_opt.point_id==1;
  opt_input.point_id = ifa_opt.point_id(unzipped);
  if zipper
    z_temp = ifa(z(unzipped), mapdata.c_array(q), mapdata.a_array(q), opt_input);
  else
    z_temp = ifa(z(unzipped), mapdata.a_array(q), opt_input);
  end
  % These points changed from being unzipped to being zipped_up
  newly_zipped = abs(imag(z_temp))>0;
  now_zipped = unzipped;
  temp = unzipped(unzipped);
  temp(not(newly_zipped)) = 0;
  now_zipped(unzipped) = temp;

  ifa_opt.point_id(now_zipped) = 0;
  z_tooth_indicator(now_zipped) = q;
  z_tooth_indicator_bins{q} = find(now_zipped);
  z(unzipped) = z_temp;
end
z(ifa_opt.point_id==1) = -z(ifa_opt.point_id==1);

% Ok, now everything is zipped up. Time to unzipper it all. 
unzip = false(size(z));
unzip(ifa_opt.point_id==1) = true; % Some points are still on R
ifa_input.point_id = ones(size(z));
ifa_input.point_id(unzip) = 2;
for q = 1:(mapdata.N_teeth)
  % Add on points that need to be unzippered
  unzip(z_tooth_indicator_bins{q}) = true;
  opt_input.point_id = ifa_input.point_id(unzip);
  if zipper
    [interior,exterior] = fa(z(unzip), mapdata.c_array(q), mapdata.a_array(q), opt_input);
  else
    [interior,exterior] = fa(z(unzip), mapdata.a_array(q), opt_input);
  end

  % Distribute interior/exterior points to where they're supposed to go
  z(unzip & zext_id) = interior(zext_id(unzip));
  z(unzip & zint_id) = exterior(zint_id(unzip));
  ifa_input.point_id(unzip) = 2; % all the unzipped points are on \mathbb{R}
end

% If zipper, do one more geodesic unzipping
if zipper
  unzip(z_tooth_indicator_bins{mapdata.N_teeth+1}) = true;
  opt_input.point_id = ifa_input.point_id(unzip);
  [interior,exterior] = fa_geo(z(unzip), mapdata.c_array(end), opt_input);

  % Distribute interior/exterior points to where they're supposed to go
  z(unzip & zext_id) = interior(zext_id(unzip));
  z(unzip & zint_id) = exterior(zint_id(unzip));
  ifa_input.point_id(unzip) = 2; % all the unzipped points are on \mathbb{R}
end

% Ok, now reapply terminal maps
% moebius:
m = [1, 0;...
     -mapdata.a_array(end), 1];
z = moebius(z,m);

z = -sign(mapdata.winding_number)*real(z.^2); % everything should be real now

% Now interior is on exterior and vice versa
z(zint_id) = moebius(z(zint_id), mapdata.m_out);
z(zext_id) = moebius(z(zext_id), mapdata.m_in);

% Ok, map to unit circle
m = [-1, i;...
      1, i];
z = moebius(z,m);
