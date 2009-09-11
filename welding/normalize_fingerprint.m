function[tin,tout] = normalize_fingerprint(tin,tout,varargin)
% [theta_int,theta_ext] = normalize_fingerprint(theta_int,theta_ext,{norm_point=1})
%
%     A very simple function: all zipper-type conformal maps are normalized so
%     that (theta_int, theta_ext) = (0,0) is a point on the fingerprint. This
%     function merely 'lines up' the values of theta_int, theta_ext so that
%     theta_int(norm_point) = 0 and theta_ext(norm_point) = 0. In the process,
%     it also unwraps the data to make it monotonic.
%
%     Basically this function just saves typing.

global handles;
opt = handles.common.input_schema({'norm_point'}, {1}, [], varargin{:});

tin = unwrap(tin);
tout = unwrap(tout);

tin = tin - tin(opt.norm_point);
tout = tout - tout(opt.norm_point);
