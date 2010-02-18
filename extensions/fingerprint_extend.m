function[h] = fingerprint_extend(stuff, x)
% fingerprint_extend -- Evaluates the fingerprint extended quasiperiodically
%
% h = fingerprint_extend(x, stuff)
%
%     Given the fingerprint piecewise polynomial construct in stuff (see
%     weld_primitive_setup), this function evaluates the quasiperiodic
%     extension of the fingerprint defined on the whole real line.

persistent ejpoly sort_histc
if isempty(ejpoly)
  from speclab.orthopoly1d.jacobi.eval import eval_jacobi_poly as ejpoly
  from labtools import sort_histc
end

xsize = size(x);
x = x(:);

shifts = floor((x/(2*pi)));
x = x - shifts*2*pi;
h = zeros(size(x));

N = length(stuff.grid.local_nodes);
% Now take the global locations and map them to [-1, 1] local locations

edges = [stuff.grid.cells(:,1); stuff.grid.cells(end,2)];
[x, bin, count, slice_matrix, sort_order] = sort_histc(x, edges);

% If something is in count(end), we've just got to place it into the
% count(end-1) slot
slice_matrix(end-1,2) = slice_matrix(end,2);
slice_matrix(end,:) = [];
bin(bin==length(count)) = length(count)-1;
count(end-1) = count(end-1) + count(end);
count(end) = [];

x = 2*(x - stuff.grid.cell_centers(bin))./stuff.grid.cell_lengths(bin);
V = ejpoly(x, 0:N-1, 'alpha', 0, 'beta', 0);

for q = 1:length(count)
  inds = slice_matrix(q,1):slice_matrix(q,2);
  if any(inds)
    h(inds) = V(inds,:)*stuff.poly_modes(:,q);
  end
end

h(sort_order) = h;
h = reshape(h + 2*pi*shifts, xsize);
