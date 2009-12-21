function[D] = delta_angle_matrix(N, varargin)
% delta_angle_matrix -- sparse matrix representation of delta operator
%
% D = delta_angle_matrix(N, {k=0, scale=dyadic})
%
%     Returns an N x N sparse matrix that performs the delta angle computation
%     on a size-N column vector z. (See delta_angles)

persistent strict_inputs indexing
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.curves import indexing
end

inputs = {'scale', 'k'};
defaults = {'dyadic', 0};

opt = strict_inputs(inputs, defaults, [], varargin{:});

ind = indexing(opt.scale);

ik = ind(opt.k);

D = -2*speye(N) + circshift(speye(N), [0,ik]) + circshift(speye(N), [0,-ik]);
D = D/(2*ik);
