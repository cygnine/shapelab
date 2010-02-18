function[stuff] = weld_primitive_setup(theta_int, theta_ext, varargin)
% weld_primitive_setup -- precomputations for primitive calculations
%
% stuff = weld_primitive_setup(theta_int, theta_ext, {local_nodes=[]})
%
%     Given the pointwise (theta_int, theta_ext) representation for the
%     fingerprint, this function interpolates those values on each local_nodes
%     stencil (where the cells are delineated by the theta_int values and
%     local_nodes is a vectors of values in [-1, 1]) and uses the resulting
%     polynomial approximation to compute various quantities for calculating the
%     primitive.

persistent strict_inputs make_pwpoly ejpoly spdiag
persistent jacobi_glq pw_primitive dcoeffs sepmat
if isempty(make_pwpoly)
  from labtools import strict_inputs spdiag
  from shapelab.zipper import make_pwpoly
  from speclab.orthopoly1d.jacobi.quad import gauss_lobatto_quadrature as jacobi_glq
  from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly as ejpoly
  from speclab.orthopoly1d.jacobi.operators import primitive_pointwise as pw_primitive
  from speclab.orthopoly1d.jacobi.coefficients import derivative as dcoeffs
  from speclab.orthopoly1d.jacobi.connection import ...
    integer_separation_connection_matrix as sepmat
end

opt = strict_inputs({'local_nodes'}, {jacobi_glq(9)}, [], varargin{:});

[p, grid] = make_pwpoly(theta_int, theta_ext, 'local_nodes', opt.local_nodes);
grid.cell_lengths = diff(grid.cells, 1, 2);
grid.cell_centers = mean(grid.cells, 2);

jopt.alpha = 0; jopt.beta = 0;

% Compute the global primitive and cellwise primitives
stuff.cell_primitives = pw_primitive(grid.local_nodes, 1, jopt)*p;
stuff.cell_primitives = (grid.cell_lengths/2)'.*stuff.cell_primitives;
stuff.global_primitive = sum(stuff.cell_primitives);
stuff.cumulative_primitives = [0, cumsum(stuff.cell_primitives(1:end-1))];

% Now compute the primitive polynomial for interval evaluation
N = length(grid.local_nodes);
V = ejpoly(grid.local_nodes, 0:(N-1), jopt);
% promote modal coefficients to (alpha+1, beta+1) class
C = sepmat(N, 0, 0, 1,1);
D = spalloc(N+1,N,nnz(C));

% normalize coefficients to correspond to (alpha,beta) multipliers
etas = dcoeffs(1:N, 0, 0);
C = spdiag(1./etas)*C;
D(2:end,:) = C;

% Find condition to fix left-endpoint
temp = speye(N+1);
temp(1,1) = 0;
temp(1,2:end) = -ejpoly(-1, 1:N, jopt)/ejpoly(-1,0,jopt);

% Yay linear algebra
integrator = temp*D*inv(V);

% This is a polynomial defined on LOCAL cell evaluations
stuff.poly = p;
stuff.poly_modes = inv(V)*p;
stuff.primitive_poly = integrator*p;

stuff.grid = grid;  % will need to project global to local
