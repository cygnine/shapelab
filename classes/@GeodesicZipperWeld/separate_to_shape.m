function[z_int, z_ext] = separate_to_shape(self, theta_int, theta_ext);
% separate_to_shape -- Runs fingerprint samples back to shape
%
% z = separate_to_shape(self, theta_int, theta_ext)
%
%     Given samples on the interior (theta_int) and/or exterior (theta_ext) of
%     the fingerprint (that is domain + range, respectively), this function runs
%     these samples backwards through the welding process onto the shape.

theta_int = theta_int(:);
theta_ext = theta_ext(:);

z_int = exp(i*theta_int);
z_ext = exp(i*theta_ext);

[z_int, z_ext] = self.inverse_moebius_alignment(z_int, z_ext);

N_interior = length(z_int);
z = [z_int; z_ext];
interior_ind = false(size(z));
interior_ind(1:N_interior) = true;
exterior_ind = not(interior_ind);

z = self.inverse_terminal_map(z, false(size(z)), interior_ind, exterior_ind);

for q = self.N_teeth:-1:1
  z = self.slider('zipup', q, z);
end

z = self.inverse_initial_map(z);

z_int = reshape(z(1:N_interior), size(theta_int));
z_ext = reshape(z(N_interior:end), size(theta_ext));
