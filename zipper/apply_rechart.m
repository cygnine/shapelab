function[theta_int, theta_ext] = apply_rechart(theta_int, theta_ext, M)
% apply_rechart -- Applies a fingerprint rechart
%
% [theta_int, theta_ext] = apply_rechart(theta_int, theta_ext, M)
%
%     Using the given 2x2 PSL_2(R) matrix M, this function applies the relevant
%     recharting to the fingerprint defined by the points (theta_int,
%     theta_ext).
%     Right now this is only done over [0, 2*pi]. I'm in no mood currently to do
%     bean counting.

persistent moebius wrap
if isempty(moebius)
  from shapelab.common import moebius
  from labtools import interval_wrap as wrap
end
interval = [0, 2*pi];

theta_int = wrap(angle(moebius(exp(i*theta_int), M)), interval);
