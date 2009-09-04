function[varargout] = interpolate_fingerprint(mapdata,varargin)
% [{theta_int_mapped, theta_ext_mapped}] = interpolate_fingerprint(mapdata,{theta_int=[], theta_ext=[]})
%
% [theta_int_mapped] = interpolate_fingerprint(mapdata, 'theta_int', theta_int)
% [theta_ext_mapped] = interpolate_fingerprint(mapdata, 'theta_ext', theta_ext)
% [theta_int_mapped, theta_ext_mapped] = interpolate_fingerprint(mapdata, 'theta_int', theta_int, 'theta_ext', theta_ext)
%
%     This function 'interpolates' a fingerprint. The first calling syntax
%     interpolates the values of theta_ext given theta_int. The second syntax
%     does vice versa. The third does both simultaneously.
%
%     This function doesn't do any 'unwrapping' of angle values; it just spits
%     out whatever Matlab's angle function gives it. 

global handles;
welding = handles.shapelab.welding;
inputs = {'theta_int', 'theta_ext'};
defaults = {[], []};
opt = handles.common.InputSchema(inputs, defaults, [], varargin{:});

s_int = size(opt.theta_int);
s_ext = size(opt.theta_ext);

N_int = prod(s_int);
N_ext = prod(s_ext);

opt.theta_int = opt.theta_int(:);
opt.theta_ext = opt.theta_ext(:);
z = [exp(i*opt.theta_int); exp(i*opt.theta_ext)];

point_id = ones([N_int + N_ext, 1]);
if N_ext>0
  point_id(N_int+1:end) = 2;
end

out = welding.switch_zipper_side(z, mapdata,'point_id', point_id);

if N_int==0;
  varargout{1} = reshape(angle(out(N_int+1:end)), s_ext);
elseif N_ext==0;
  varargout{1} = reshape(angle(out(1:N_int)), s_int);
else
  varargout{1} = reshape(angle(out(1:N_int)), s_int);
  varargout{2} = reshape(angle(out(N_int+1:end)), s_ext);
end
