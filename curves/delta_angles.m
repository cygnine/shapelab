function[delta] = delta_angles(z,varargin)
% delta_angles -- computation of delta angles from complex samples
%
% delta = delta_angles(z, {k=0, scale=dyadic})
%
%     Computes the delta angles at scale k, defined as 
%
%        delta_k(n) = z(n+ind(k)) - 2*z(n) + z(n-ind(k)),
%
%     for all n. ind(k) is an indexing function defined depending on the
%     optional input 'scale':
%
%        scale = 'dyadic'            Dyadic scaling: ind(k) = 2^k
%        scale = 'linear'            Linear scaling: ind(k) = k+1
%
%     The range of the angle is the same as Matlab's `angle' function: [-pi,pi].
%     The delta angles are linear approximations to the beta angles (see
%     beta_angles). Note that the normalization here assumes that |z(n+1) -
%     z(n)| ~ 1. 

persistent strict_inputs indexing
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.curves import indexing
end

inputs = {'scale', 'k'};
defaults = {'dyadic', 0};

opt = strict_inputs(inputs, defaults, [], varargin{:});

zsize = size(z);
z = z(:);

ind = indexing(opt.scale);
ik = ind(opt.k);

zp = circshift(z, -ik);
zm = circshift(z, ik);

delta = (zp+zm-2*z)/(2*ik);

delta = reshape(delta, zsize);
