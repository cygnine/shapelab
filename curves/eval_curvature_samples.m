function[r] = eval_curvature_samples(k,s)
% eval_curvature_samples -- Evaluates a plane curve by specifying curvature
%
% r = eval_curvature_samples(k,s)
%
%     With a curve specified by its polynomial curvature coefficients k, this
%     function evaluates the curve by sampling the curve at the arc-length
%     points s. These points should be monotonically increasing.
%
%     The returned quantity r is a complex-valued vector of the same size as s.

persistent trig_integral
if isempty(trig_integral)
  from shapelab.curves import trig_integral
end

ssize = size(s);

s = s(:);

s = [[0; s(1:(end-1))] s];

[cint,sint] = trig_integral(k,s);

r = cumsum(cint + i*sint);
r = reshape(r, ssize);
