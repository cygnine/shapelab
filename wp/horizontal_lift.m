function[vz,m,w] = horizontal_lift(nodes, v, z, d)
% horizontal_lift -- Evaluates the WP horizontal lift of velocity fields
%
% [vz,m,w] = horizontal_lift(nodes, v, z, [[d=0]])
%
%     If nodes and v are vectors of the same size, this function computes the WP
%     horizontal lift coefficients m and w that correspond to the velocity field
%     evaluations (nodes, v). The horitzontal lift is evaluated at the locations
%     z, and fz is of the same size as z. The d'th derivative is evaluated, and
%     so far this function is valid for d=1,2. If N = length(v), then N =
%     length(m) and 3 = length(w).
%
%     If nodes and v are arrays of the same size, then the input z must be an
%     array with the same number of columns. Let [N,M] = size(nodes). Then this
%     function assumes that there are M collections of length-N nodes and data:
%     (nodes(:,m), data(:,m)). The horizontal lift of all of these is computed
%     by assigning equal weight to each lift and finding the kernel coefficients
%     w that simultaneously minimize the sum of WP norms of the lifts. Then
%     [N,M] = size(m) and 3 = length(w). Then for m = 1, 2, ..., M, the m'th
%     horizontal lift is evaluated at z(:,m) and the result is stored in
%     vz(:,m). Again, the d'th derivative is computed.
%
%     The input z can be empty in which case vz is empty and m and w contain the
%     coefficients.
%
%     The workhorse for this routine is the gsvd.

persistent greens_matrix kernel_basis greens_function
if isempty(greens_matrix)
  from shapelab.wp import greens_matrix kernel_basis greens_function
end

[N,M] = size(nodes);

% Input validation
if nargin < 4
  d = 0;
  if nargin < 3
    z = [];
  end
end
if any([N,M] ~= size(v));
  error('Inputs ''nodes'' and ''v'' must be of the same size');
end
if (M ~= size(z,2)) & (size(z,2) > 0)
  error('Input ''z'' must have the same number of columns as ''nodes''');
end

% These can be huge matrices:
% Possible TODO: sparse matrix with svds --> gsvds?
G = zeros([N*M N*M+3]);
L = zeros([N*M N*M+3]);

for m = 1:M
  i1 = 1 + (m-1)*N;
  i2 = m*N;
  G(i1:i2,i1:i2) = greens_matrix(nodes(:,m));
  L(i1:i2,i1:i2) = chol(G(i1:i2,i1:i2));

  G(i1:i2,end-2:end) = kernel_basis(nodes(:,m));
end

[U,V,X,C,S] = gsvd(G, L);

% The old-school way
%pinvC = spdiags(1./diag(C,3), -3, N*M+3, N*M);
%H = inv(X')*pinvC*U';
%mw = H*v(:);

% All the cool kids are doing it this equivalent way;
diagC = 1./diag(C,3);
mw = X'\[zeros([3 1]); diagC.*(U'*v(:))];

mom = reshape(mw(1:M*N), [N M]);
w = mw(M*N+1:end);

% Now evaluate the lift
if isempty(z)
  vz = [];
else
  vz = zeros(size(z));
  for m = 1:M
    [x,y] = meshgrid(z(:,m), nodes(:,m));
    vz(:,m) = greens_function((x-y)', d)*mom(:,m) + ...
              kernel_basis(z(:,m), d)*w;
  end
end

m = mom;
