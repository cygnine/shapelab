function[c] = wp_minimial_ensemble(theta, v, varargin)
% [c] = wp_minimal_ensemble(theta, v, {Ginvs = [], Us = []})
% 
%     This function computes the element of the kernel of the WP operator that
%     *simultaneously* minimizes the inv(G)-norm for the data
%
%      theta(:,1)   v(:,1)
%      theta(:,2)   v(:,2)
%      ...................
%      theta(:,end) v(:,end)
%
%     In other words, the minimal value of the metric is attained at \sum_k
%     w(:,k)'*inv(G_k)*w(:,k), where w(:,k) = v(:,k) - U_k*c, and U is the
%     coordinate map between a basis of the WP operator and point-evaluations on
%     theta(:,k).
%
%     If the cell-array Ginvs contains all the inverted G-matrices, this saves
%     time.

persistent strict_inputs kernel_basis greens_function
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.wp import greens_function kernel_basis
end

opt = strict_inputs({'Ginvs', 'Us'}, {[], []}, [], varargin{:});
K = size(theta, 2);

if isempty(opt.Ginvs)
  opt.Ginvs = cell([K 1]);
  for q = 1:K
    [x,y] = mesh_grid(theta(:,q), theta(:,q));
    opt.Ginvs{q} = inv(greens_function(x-y));
  end
end

if isempty(opt.Us)
  opt.Us = cell([K 1]);
  for q = 1:K
    opt.Us{q} = kernel_basis(theta(:,q));
  end
end

Q = zeros(3);
rhs = zeros([3 1]);
for q = 1:K
  temp = opt.Us{q}'*opt.Ginvs{q};
  Q = Q + temp*opt.Us{q};
  rhs = rhs + temp*v(:,q);
end

c = inv(Q)*rhs;
