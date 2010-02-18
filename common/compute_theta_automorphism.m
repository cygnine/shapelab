function[H] = compute_theta_automorphism(theta_in, theta_out)
% compute_theta_automorphism -- Returns a member of PSL_2(R)
%
% H = compute_theta_automorphism(theta_in, theta_out)
%
%     Returns a member of PSL_2(R) that maps h(theta_in) to h(theta_out), where
%     h(theta) = tan(theta/2). The implementation of this map is in
%     theta_autormorphism.

persistent find_map
if isempty(find_map)
  from shapelab.common.moebius_maps import psl2_as_moebius as find_map
end

z1 = tan(theta_in/2);
z2 = tan(theta_out/2);

H = find_map(z1, z2);
