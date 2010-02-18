function[p, grid] = make_pwpoly(theta_int, theta_ext, varargin)
% make_pwpoly -- Creates a piecewise polynomial from a fingerprint
%
% [p, grid] = make_pwpoly(theta_int, theta_ext, {local_nodes = []})
%
%     Given K points (theta_int, theta_ext) that define a fingerprint, this
%     function defines an element between two theta_int values and uses
%     the local_nodes (N values in [-1, 1]) for evaluation at those locations.
%     This is a classic finite-element setup. The nodal evaluations are
%     returned in the N x K matrix p.
%
%     The fingerprint is interpolated using the geodesic version of the
%     Zipper algorithm. The default local_nodes are 4
%     Legendre-Gauss-Lobatto points. If you want to form a function
%     f(theta_ext) ----> theta_int instead, just reverse the inputs.
%
%     TODO: fix theta_int(1) \neq 0 and theta_ext(1) \neq 0 bugs

persistent interp_fprint ejpoly strict_inputs jacobi_glq repnodes
persistent interval_wrap
if isempty(strict_inputs)
  from labtools import strict_inputs interval_wrap
  from shapelab.zipper import interpolate_fingerprint as interp_fprint
%  from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly as ejpoly
  from speclab.orthopoly1d.jacobi.quad import gauss_lobatto_quadrature as jacobi_glq
  from dg.meshes import replicate_local_nodes as repnodes
end

opt = strict_inputs({'local_nodes'}, {jacobi_glq(5)}, [], varargin{:});
theta_int = interval_wrap(theta_int(:), [0, 2*pi]);
if isempty(opt.local_nodes)
  opt.local_nodes = jacobi_glq(5);
end

N = length(opt.local_nodes);
%V = ejpoly(opt.local_nodes, 0:(N-1));
cells = [theta_int, [theta_int(2:end); theta_int(1) + 2*pi]];
p = repnodes(opt.local_nodes, cells);

temp = interp_fprint(theta_int, theta_ext, 'theta_int', p);
%temp = interp_fprint(theta_int, theta_ext, 'theta_int', p(:,1));
% Stupid machine eps crap:
if abs(p(1,1))<1e-14;
  temp(1,1) = 0;
end
if abs(p(end,end)-2*pi)<1e-14
  temp(end,end) = 2*pi;
end
p = temp;
%p = inv(V)*temp;

grid.cells = cells;
grid.local_nodes = opt.local_nodes;
