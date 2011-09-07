function[w] = wp_minimal_interpolation(theta, v, phi)
% wp_minimal_interpolation -- Interpolates the data with a WP-minimal lift
%
% w = wp_minimal_interpolation(theta, v, phi)
%
%     Given the data (theta, v), where each entry of theta lies in [0, 2*pi),
%     this function interpolates the function defined by (theta,v) at the
%     locations phi on [0, 2*pi). The constructed function is the WP-minimal
%     horizontal lift -- see wp_minimal_metric for details.

persistent wp_minimal_metric kernel_basis greens_function
if isempty(wp_minimal_metric)
  from shapelab.wp import wp_minimal_metric kernel_basis greens_function
end

phisize = size(phi);
phi = phi(:);
v = v(:);

[H, K] = wp_minimal_metric(theta);
[temp1, temp2] = meshgrid(theta, phi);
G = greens_function(temp1-temp2);
U = kernel_basis(phi);

w = G*(H*v) + U*(K*v);

w = reshape(w, phisize);
