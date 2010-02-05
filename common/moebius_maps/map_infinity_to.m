function[H] = map_infinity_to(image)
% map_infinity_to -- Generates a Moebius map based on the image of infinity
%
% H = map_infinity_to(image)
%
%     Given a scalar real-valued point image, this returns the matrix
%     representation of the Moebius map that is the pure rotation of S^2 about
%     the plane generated the geodesic passing through the image point and
%     infinity. The point infinity is mapped to the image.

persistent conjugate_points
if isempty(conjugate_points)
  from shapelab.common.moebius_maps import conjugate_points
end

if image==0
  H = [0 1; 1 0]; 
  return
end

% First find the conjugate points to the geodesic
[mu, nu] = conjugate_points(image, Inf);

A = [ones([2, 1]), -[mu; nu]];
coeffs = inv(A)*[mu^2-image*mu; nu^2-image*nu];

H = [image, coeffs(1); 1, coeffs(2)];
