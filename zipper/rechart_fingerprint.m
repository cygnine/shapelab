function[vertices_int, vertices_ext] = rechart_fingerprint(vertices_int, vertices_ext, theta_int, theta_ext);
% rechart_fingerprint -- Recharts a fingerprint
%
% [vertices_int, vertices_ext] = rechart_fingerprint(vertices_int, vertices_ext, theta_int, theta_ext)
%
%     Given point values (vertices_int, vertices_ext) that define a fingerprint,
%     this function applies an interior Moebius map to the fingerprint so that
%     (theta_int(i), theta_ext(i)) is a point on the fingerprint for each i =
%     1,2,3.
%
%     Right now this is only done over [0, 2*pi]. I'm in no mood currently to do
%     bean counting.

persistent interp_fprint find_moebius moebius wrap
if isempty(interp_fprint)
  from shapelab.zipper import interpolate_fingerprint as interp_fprint
  from shapelab.common.moebius_maps import specify_points as find_moebius
  from shapelab.common import moebius
  from labtools import interval_wrap as wrap
end

% First send theta_ext through the fingerprint, and use a Moebius map to take
% whatever comes out to the locations theta_int
[garbage, current_theta_int] = interp_fprint(vertices_int, vertices_ext, 'theta_ext', theta_ext);

z1 = exp(i*current_theta_int);
z2 = exp(i*theta_int);

M = find_moebius(z1, z2);

interval = [0, 2*pi];

vertices_int = wrap(angle(moebius(exp(i*vertices_int), M)), interval);
