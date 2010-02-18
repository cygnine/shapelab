function[psi] = weld_primitive_evaluate(theta, stuff)
% weld_primitive_evaluate -- Evaluates the primitive of a welding map
%
% psi = weld_primitive_evaluate(theta, stuff)
%
%     This function uses the information about the primitive in the struct stuff
%     to evaluate the primitive at the locations theta (which must all be
%     in [0, 2*pi]). No checking is performed to make sure that theta lies in
%     this interval. The struct stuff comes from the output of
%     weld_primitive_setup.

tsize = size(theta);
theta = theta(:);

persistent ejpoly sort_histc
if isempty(ejpoly)
  from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly as ejpoly
  from labtools import sort_histc
end

% First we use histc to locate where the inputs are
edges = [stuff.grid.cells(:,1); stuff.grid.cells(end,2)];
[theta, bin, count, slice_matrix, sort_order] = sort_histc(theta, edges);

% If something is in count(end), we've just got to place it into the
% count(end-1) slot
slice_matrix(end-1,2) = slice_matrix(end,2);
slice_matrix(end,:) = [];
bin(bin==length(count)) = length(count)-1;
count(end-1) = count(end-1) + count(end);
count(end) = [];

% length(count) is the number of bins

% Now take the global locations and map them to [-1, 1] local locations
theta = 2*(theta - stuff.grid.cell_centers(bin))./stuff.grid.cell_lengths(bin);

% Sum up precomputed primitives:
psi = stuff.cumulative_primitives(bin).';

N = length(stuff.grid.local_nodes);
% Generate a global Vandermonde matrix
V = ejpoly(theta, 0:N, 'alpha', 0, 'beta', 0);

for q = 1:length(count);
  inds = slice_matrix(q,1):slice_matrix(q,2);
  if any(inds)
    psi(inds) = psi(inds) + stuff.grid.cell_lengths(q)/2*...
                            V(inds,:)*stuff.primitive_poly(:,q);
  end
end

% Unsort
psi(sort_order) = psi;
psi = reshape(psi, tsize);
