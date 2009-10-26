% This examples takes a 4-leaf clover and does some fingerprint stuff with it.

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

% First construct points on the clover:
clover_opt.lobes = 4;  % 4-leafed
clover_opt.lobe_depth_ratio = 0.95; % Large lobes, collapsing nearly to 0
N = 100; % number of vertices
z = clover(N, clover_opt);

% FYI, the 'node_locations' optional input to clover lets you put points inside
% or outside the clover. Not really useful in this example, though.

% For no good reason, plot the clover
figure(1);
plot(z, 'b.'); hold on;
ltex(title('A 4-leafed clover'));

% Ok, now compute three different maps:
%  - a purely 'geodesic' map
%  - a purely 'slit' map
%  - a 'zipper' map + terminal geodesic map (necessary for welding)

% Both of these are optional: the pre-images of 0 and infinity in shape-space
map_opt.shape_0 = 0;
map_opt.shape_infinity = Inf;

map_opt.type = 'geodesic';
fprintf('Computing geodesic map....\n');
geo_map = welding.make_zipper_map(z,map_opt);

map_opt.type = 'slit';
fprintf('Computing slit map....\n');
slit_map = welding.make_zipper_map(z,map_opt);

map_opt.type = 'zipper_weld'; % NOTE: DO NOT use 'zipper' ... that *is* a
                              % mapping algorithm programmed, but it's not 
                              % weld-able
fprintf('Computing zipper map....\n');
zip_map = welding.make_zipper_map(z,map_opt);

% You can inspect the three map structs if you like...you'll probably be able to
% figure some stuff out...but it's not that telling. These structs basically
% completely define a map from the interior of the shape to the unit disc, and
% from the exterior of the shape to complement of the unit disc. On the unit
% disc, the function is multi-valued, which you will see causes headaches when
% telling the computer what you want it to do.

% First thing: each of the structs saves the unzipped (interior + exterior)
% locations, so the pointwise fingerprint at the vertices is already available.
% Unwrap
tempf = @(x) unwrap(x);
tin_geo = tempf(angle(geo_map.vertices_in));
tout_geo = tempf(angle(geo_map.vertices_out));

tin_slit = tempf(angle(slit_map.vertices_in));
tout_slit = tempf(angle(slit_map.vertices_out));

tin_zip = tempf(angle(zip_map.vertices_in));
tout_zip = tempf(angle(zip_map.vertices_out));

figure(2);
subplot(2,2,1);
plot(tin_geo, tout_geo, 'b.');
ltex(xlabel('$\theta_{in}$'));
ltex(ylabel('$\theta_{out}$'));
ltex(title('Geodesic fingerprint'));
subplot(2,2,2);
plot(tin_slit, tout_slit, 'b.');
ltex(xlabel('$\theta_{in}$'));
ltex(ylabel('$\theta_{out}$'));
ltex(title('Slit fingerprint'));
subplot(2,2,3);
plot(tin_zip, tout_zip, 'b.');
ltex(xlabel('$\theta_{in}$'));
ltex(ylabel('$\theta_{out}$'));
ltex(title('Zipper fingerprint'));
subplot(2,2,4)
plot(tin_geo, tout_geo, 'b.', tin_slit, tout_slit, 'r+', ...
     tin_zip, tout_zip, 'ko');
ltex(xlabel('$\theta_{in}$'));
ltex(ylabel('$\theta_{out}$'));
ltex(title('All fingerprints'));
legend('Geodesic', 'Slit', 'Zipper');

% You can inspect for yourself that all of the functions generated are
% monotonic; the numerical stability is actually quite good regarding that.

% You see from the last quadrant that comparing them is difficult because they
% are all at different heights. Wouldn't it be nice if we could interpolate them
% all to the same points and then subtract out the (hopefully constant)
% difference to really compare them? Later...right now:

% Since we know the map, we can invert it, so let's put a crapload of points
% onto the unit circle and invert the map and see exactly what shape we're
% approximating with each map. In this case, we have to tell the algorithm
% whether we're on the interior of the unit circle or the exterior; the
% map has different values for each case.
w_fine = shapelab.common.polar_linspace(1,100*N,'r0', 1, 'r1', 1);

% The optional input 'point_id' tells the inverse map where the points are. See
% that function's helpstring
point_id = ones(size(w_fine));  % interior of unit circle
z_geo_fine = zipper.evaluate_inverse_map(w_fine, geo_map,'point_id', point_id);
z_slit_fine = zipper.evaluate_inverse_map(w_fine, slit_map,'point_id', point_id);
z_zip_fine = zipper.evaluate_inverse_map(w_fine, zip_map,'point_id', point_id);

% Why not combine these with the points from the exterior map as well?
point_id = 2*ones(size(w_fine)); % exterior of unit circle
z_geo_fine = [z_geo_fine; zipper.evaluate_inverse_map(w_fine, geo_map,'point_id', point_id)];
z_slit_fine = [z_slit_fine; zipper.evaluate_inverse_map(w_fine, slit_map,'point_id', point_id)];
z_zip_fine = [z_zip_fine; zipper.evaluate_inverse_map(w_fine, zip_map,'point_id', point_id)];

% We'll also have a fine-scale plot of the original clover on top as well.
z_fine = clover(100*N,clover_opt);

figure(3);
temp = [-1.05 1.05 -1.05 1.05];
subplot(2,2,1);
title('The geodesic approximation')
plot(z_geo_fine, 'r.'); hold on; plot(z_fine, 'b'); set(plot(z, 'k.'), 'markersize', 10);
axis(temp);
subplot(2,2,2);
title('The slit approximation')
plot(z_slit_fine, 'r.'); hold on; plot(z_fine, 'b'); set(plot(z, 'k.'), 'markersize', 10);
axis(temp);
subplot(2,2,3);
title('The zipper approximation')
plot(z_zip_fine, 'r.'); hold on; plot(z_fine, 'b'); set(plot(z, 'k.'), 'markersize', 10);
axis(temp);
subplot(2,2,4);
plot(0,0, 'r.'); hold on; plot(0,0,'b'); hold on; plot(0,0,'k.');
axis([1,2,1,2]);
legend('Mapped inverse values', 'Original shape', 'Original data points');

% Looking at this previous graph, you can tell that the approximations indeed
% are different. It's most apparent near the origin. The red dots represent the
% real shape that the map defines. Note that because of this, it's numerically
% impossible to have a "forward" mapping from the shape onto the unit circle:
% we'd have to put points exactly on that weird shape, which we don't know.

% Just for fun, however, we can also invert the map from the interior and
% exterior of the unit circle. A forward map *is* possible from the
% exterior/interior of the shape to the complement of the unit circle, but I
% haven't implemented it due to the fact that I don't see the pressing utility
% of it.

% The following just creates a refined mesh of the unit circle.
w_int_fine = [];
w_ext_fine = [];
N_stages = 6;  % number of refinements
refine_int_opt.r0 = 0;
refine_int_opt.r1 = 0.75;

refine_ext_opt.r1 = 3;
refine_ext_opt.r0 = 1.25;
for q = 1:N_stages
  w_int_fine = [w_int_fine; shapelab.common.polar_linspace(20, 2^(q-1)*N/4, refine_int_opt)];
  refine_int_opt.r0 = refine_int_opt.r1;
  refine_int_opt.r1 = 1 - (0.25)^(q+1);

  w_ext_fine = [w_ext_fine; shapelab.common.polar_linspace(20, 2^(q-1)*N, refine_ext_opt)];
  refine_ext_opt.r1 = refine_ext_opt.r0;
  refine_ext_opt.r0 = 1 + (0.25)^(q+1);
end

fprintf(['Computing interior/exterior inverse evaluations (I chose a lot of points, \n'...
         'so this takes some time....\n']);

point_id = zeros(size(w_int_fine));  % anywhere not on the unit circle (this is the
                                     % default behavior)
z_geo_int_fine = zipper.evaluate_inverse_map(w_int_fine, geo_map,'point_id', point_id);
z_slit_int_fine = zipper.evaluate_inverse_map(w_int_fine, slit_map,'point_id', point_id);
z_zip_int_fine = zipper.evaluate_inverse_map(w_int_fine, zip_map,'point_id', point_id);

point_id = zeros(size(w_ext_fine));  % anywhere not on the unit circle (this is the
                                     % default behavior)
z_geo_ext_fine = zipper.evaluate_inverse_map(w_ext_fine, geo_map,'point_id', point_id);
z_slit_ext_fine = zipper.evaluate_inverse_map(w_ext_fine, slit_map,'point_id',point_id);
z_zip_ext_fine = zipper.evaluate_inverse_map(w_ext_fine, zip_map,'point_id', point_id);

% No reason, just because I like purty pictures
figure(4);
temp = [-2, 2, -2, 2];
subplot(2,2,1);
plot(z_geo_fine, 'r.'); hold on; 
set(plot(z_geo_int_fine, 'b.'), 'markersize', 2);
set(plot(z_geo_ext_fine, 'k.'), 'markersize', 2);
axis(temp);
ltex(title('Geodesic approximation'));
subplot(2,2,2);
plot(z_slit_fine, 'r.'); hold on; 
set(plot(z_slit_int_fine, 'b.'), 'markersize', 2);
set(plot(z_slit_ext_fine, 'k.'), 'markersize', 2);
ltex(title('Slit approximation'));
axis(temp);
subplot(2,2,3);
plot(z_zip_fine, 'r.'); hold on; 
set(plot(z_zip_int_fine, 'b.'), 'markersize', 2); 
set(plot(z_zip_ext_fine, 'k.'), 'markersize', 2);
ltex(title('Zipper approximation'));
axis(temp);

% Ok, to the punch line. Let's "interpolate" the fingerprints are particular
% locations. We can do this in the "Riemann" sense (stipulating tin), or the
% "Lebesgue" sense (stipulating tout); or a combination of both if you want. 

fprintf('Interpolating fingerprint values....\n');
% Let's say I want tout as a function of tin for equispaced tin:
tin = linspace(0, 2*pi, 1000);

tin_geo_image = welding.interpolate_fingerprint(geo_map, 'theta_int', tin);
tin_slit_image = welding.interpolate_fingerprint(slit_map, 'theta_int', tin);
tin_zip_image = welding.interpolate_fingerprint(zip_map, 'theta_int', tin);

% Or the "Lebesgue" way...specify tout:
tout = linspace(0, 2*pi,1000);
tout_geo_image = welding.interpolate_fingerprint(geo_map, 'theta_ext', tout);
tout_slit_image = welding.interpolate_fingerprint(slit_map, 'theta_ext', tout);
tout_zip_image = welding.interpolate_fingerprint(zip_map, 'theta_ext', tout);

% Now for normalization of the fingerprints so that they line up. *All* maps are
% normalized so that (tin, tout) = (0,0) corresponds to the first vertex of the
% shape, so all these points line up. In all the inputs above (tin or tout), the
% first entry was 0, so we can normalize all these points back down to 0 to line
% everything up. There's a function that does this automatically. 
[tin, tin_geo_image] = fprint_norm(tin,tin_geo_image);
[tin, tin_slit_image] = fprint_norm(tin,tin_slit_image);
[tin, tin_zip_image] = fprint_norm(tin,tin_zip_image);

[tout_geo_image, tout] = fprint_norm(tout_geo_image,tout);
[tout_slit_image, tout] = fprint_norm(tout_slit_image,tout);
[tout_zip_image, tout] = fprint_norm(tout_zip_image,tout);

% Let's see that indeed we have interpolated the functions:
figure(5);
subplot(2,2,1);
plot(tin, tin_geo_image, 'r.', tout_geo_image, tout, 'b.', tin_geo, tout_geo, 'k*');
ltex(xlabel('$\theta_{int}$'));
ltex(ylabel('$\theta_{ext}$'));
ltex(title('Geodesic fingerprint'));
ltex(legend('Interpolating $\theta_{ext}$', 'Interpolating $\theta_{int}$', 'Original fingerprint samples'));
subplot(2,2,2);
plot(tin, tin_slit_image, 'r.', tout_slit_image, tout, 'b.', tin_slit, tout_slit, 'k*');
ltex(xlabel('$\theta_{int}$'));
ltex(ylabel('$\theta_{ext}$'));
ltex(title('Slit fingerprint'));
subplot(2,2,3);
plot(tin, tin_zip_image, 'r.', tout_zip_image, tout, 'b.', tin_zip, tout_zip, 'k*');
ltex(xlabel('$\theta_{int}$'));
ltex(ylabel('$\theta_{ext}$'));
ltex(title('Zipper fingerprint'));

% You can zoom in and inspect the fingerprints...the interpolations look quite
% good. You can also check that they're monotonic: also very good. These are the
% exact (blah blah numerical precision) values of the fingerprint of the shapes
% from figure 3. 

% Overlapping fingerprints doesn't really work yet...need to normalize them
figure(6);
subplot(2,1,1);
plot(tin, tin_geo_image, 'b.', tin, tin_slit_image, 'k.', tin, tin_zip_image, 'r.');
ltex(xlabel('$\theta_{int}$'));
ltex(ylabel('$\theta_{ext}$'));
ltex(title('Interpolating values of $\theta_{ext}$'));
set(gca, 'FontSize', 12, 'fontweight', 'b');
legend('Geodesic', 'Slit', 'Zipper');
subplot(2,1,2);
plot(tout_geo_image, tout, 'b.', tout_slit_image, tout, 'k.', tout_zip_image, tout, 'r.');
ltex(xlabel('$\theta_{int}$'));
ltex(ylabel('$\theta_{ext}$'));
ltex(title('Interpolating values of $\theta_{int}$'));
set(gca, 'FontSize', 12, 'fontweight', 'b');
legend('Geodesic', 'Slit', 'Zipper');
