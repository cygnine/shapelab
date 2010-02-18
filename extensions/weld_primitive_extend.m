function[psi] = weld_primitive_extend(theta_int, theta_ext, x, varargin)
% weld_primitive_extend -- Evaluates the real-line extension of the weld primitive
%
% psi = weld_primitive_extend(theta_int, theta_ext, x, {local_nodes=[]})
%
%     Using the fingerprint evaluations (theta_int, theta_ext), this function
%     uses the zipper + primitive calculations of this package's
%     weld_primitive series to evaluate the primitive of the extension
%     of the fingerprint to the real line at the locations x.

persistent strict_inputs wsetup weval wextend
if isempty(weval)
  from labtools import strict_inputs
  from shapelab.extensions import weld_primitive_setup as wsetup
  from shapelab.extensions import weld_primitive_evaluate as weval
  from shapelab.extensions import weld_primitive_extend_driver as wextend
end

opt = strict_inputs({'local_nodes'}, {[]}, [], varargin{:});

stuff = wsetup(theta_int, theta_ext, 'local_nodes', opt.local_nodes);

psi = wextend(x, stuff);
