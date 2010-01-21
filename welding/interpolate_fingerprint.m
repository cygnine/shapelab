function[theta_map] = interpolate_fingerprint(mapdata,theta, varargin)
% interpolate_fingerprint -- Interpolates values of the fingerprint
%
% theta_map = interpolate_fingerprint(mapdata, theta, {domain='interior'})
%
%     This function 'interpolates' a fingerprint. If domain is 'interior', then
%     the input values theta are assumed to be "interior" values, and they are
%     mapped to the corresponding exterior values (theta_map) using the welding
%     map supplied by mapdata. The opposite occurs if domain is set to
%     'exterior'.
%
%     This function doesn't do any 'unwrapping' of angle values; it just spits
%     out whatever Matlab's angle function gives it. 

persistent strict_inputs switch_zipper_side
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.welding import switch_zipper_side
end

inputs = {'domain'};
defaults = {'interior'};
opt = input_schema(inputs, defaults, [], varargin{:});

theta_size = size(theta);

z = exp(i*theta(:));

theta_map = switch_zipper_side(z, mapdata, opt);
theta_map = reshape(angle(theta_map), theta_size);
