function[beta] = beta_angles(z,varargin)
% beta_angles -- computation of beta angles from complex samples
%
% beta = beta_angles(z, {k=0, scale=dyadic})
%
%     Computes the beta angles at scale k, defined as 
%
%        beta_k(n) = arg((z(n+ind(k)) - z(n))/(z(n) - z(n-ind(k)))),
%
%     for all n. ind(k) is an indexing function defined depending on the
%     optional input 'scale':
%
%        scale = 'dyadic'            Dyadic scaling: ind(k) = 2^k
%        scale = 'linear'            Linear scaling: ind(k) = k+1
%
%     The range of the angle is the same as Matlab's `angle' function: [-pi,pi].

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

zp = circshift(z, -ind(opt.k));
zm = circshift(z, ind(opt.k));

beta = angle((zp-z)./(z-zm));

beta = reshape(beta, zsize);
