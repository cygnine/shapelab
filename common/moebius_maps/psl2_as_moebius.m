function[H] = psl2_as_moebius(z1, z2)
% psl2_as_moebius -- Returns a member from PSL_2(\mathbb{R})
%
% H = psl2_r(z1, z2)
%
%     Returns the 2x2 matrix representation of a member from PSL_2(R). The
%     member is chosen via interpretation as a Moebius map. Both the
%     3-vectors z1 and z2 are assumed to be real-valued. The map that
%     takes z1 to z2 is returned.

persistent find_map
if isempty(find_map)
  from shapelab.common.moebius_maps import specify_points as find_map
end

if nargin<3
  orientation=0;
end

H = find_map(z1, z2);
H = H/sqrt(abs(det(H)));  % So det(H) = 1

if det(H)<0
  error('These points do not define an orientiation-preserving transformation');
end
