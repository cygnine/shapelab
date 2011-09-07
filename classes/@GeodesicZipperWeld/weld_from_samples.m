function[self] = weld_from_samples(self, theta_interior, theta_exterior)
% weld_from_samples -- Computes the welding map given a collection of samples
%
% self = weld_from_samples(self, theta_interior, theta_exterior)
%
%     Given samples (theta_interior, theta_exterior) on [0, 2*pi) of the welding
%     map, this function determines the weld.
%
%     The map ensures that self.center_location is mapped to self.center_image,
%     and that self.inf_location is mapped to self.inf_image.

% The strategy here is as follows: we first backtrack through the zipper map
% without any terminal moebius alignments. Then we push self.inf_location and
% self.center_location forward through the map to determine what the terminal
% moebius maps should be. 
%
% Finally, this allows us to define the weld.

% Let's do some argument checking
theta_interior = theta_interior(:);
theta_exterior = theta_exterior(:);
%[theta_interior, inds] = sort(theta_interior);
%theta_exterior = theta_exterior(inds);

[sorted_interior, inds] = sort(theta_interior);
sorted_exterior = theta_exterior(inds);

% I don't think this check is necessary -- we do exp(i*.) anyway, and then do
% appropriate rotations
%if any(diff(sorted_exterior)<= 0) | any(diff(sorted_interior)<=0)
%  error('The input samples are not strictly monotone');
%end

self.interior_vertices = theta_interior;
self.exterior_vertices = theta_exterior;

fd_at_inf = [1e4 1e4+1].';

% The only initial moebius alignment necessary then is a rotation
z_int = [exp(i*theta_interior(:)); self.center_image];
z_ext = [exp(i*theta_exterior(:)); fd_at_inf; self.inf_image];
self.N = length(theta_interior);
self.N_teeth = self.N - 2;

self.interior_disc_vertices = z_int(1:self.N);
self.exterior_disc_vertices = z_ext(1:self.N);

[z_int, z_ext, self] = self.calculate_inverse_moebius_alignment(z_int, z_ext);

% Now z_int(1) and z_ext(1) are at Inf, and z_int(end), z_ext(end) are at 0. The
% first part of the terminal map takes a sqrt; then we need to find the
% moebius Map that takes the locations of z_ext(end-1) and z_int(end-1) to
% symmetric locations around 0. 
[z_int, z_ext, self] = self.calculate_inverse_terminal_map(z_int, z_ext);

for q = self.N_teeth:-1:1
  %[z] = self.slider('zipup', q, z, interior, slit_interior, slit_exterior);
  [z_int, z_ext, self] = self.zipup(q, z_int, z_ext);
end

[z_int, z_ext, self] = self.calculate_inverse_initial_map(z_int, z_ext);
