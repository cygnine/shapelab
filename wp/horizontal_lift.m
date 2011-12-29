function[vz,m,w,v2p,p2v] = horizontal_lift(nodes, v, z, d, ht)
% horizontal_lift -- Evaluates the WP horizontal lift of velocity fields
%
% [vz,m,w,v2p,p2v] = horizontal_lift(nodes, v, z, [[d=0, ht=[] ]])
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
%     The input 'ht' determines weights to be put on each ensemble member.
%     Therefore, length(ht)==M.
%
%     The workhorse for this routine is the gsvd. The last two outputs are
%     function handles that allow one to travel between momentum and velocity: 
%     v = v2p(p) and p = p2v(v).

persistent greens_matrix kernel_basis greens_function
if isempty(greens_matrix)
  from shapelab.wp import greens_matrix kernel_basis greens_function
end

[N,M] = size(nodes);

% Input validation
if nargin < 5
  ht = ones([M 1]);
  if nargin < 4
    d = 0;
    if nargin < 3
      z = [];
    end
  end
end
if any([N,M] ~= size(v));
  error('Inputs ''nodes'' and ''v'' must be of the same size');
end
if (M ~= size(z,2)) & (size(z,2) > 0)
  error('Input ''z'' must have the same number of columns as ''nodes''');
end
if M ~= length(ht)
  error('Input ''ht'' must have length equal to number of columns of ''nodes''');
end

% These can be huge matrices:
% Possible TODO: sparse matrix with svds --> gsvds?
G = zeros([N*M N*M+3]);
L = zeros([N*M N*M+3]);

for m = 1:M
  i1 = 1 + (m-1)*N;
  i2 = m*N;
  G(i1:i2,i1:i2) = greens_matrix(nodes(:,m));
  [vtemp,dtemp] = eig(G(i1:i2,i1:i2));
  dtempd = diag(dtemp);
  dtempd(dtempd<0) = 0;
  dtemp = sqrt(diag(dtempd));
  L(i1:i2,i1:i2) = sqrt(ht(m))*dtemp*vtemp';
  %L(i1:i2,i1:i2) = sqrt(ht(m))*chol(G(i1:i2,i1:i2));

  G(i1:i2,end-2:end) = kernel_basis(nodes(:,m));
end

[U,V,X,C,S] = gsvd(G, L);

% The old-school way
%pinvC = spdiags(1./diag(C,3), -3, N*M+3, N*M);
%H = inv(X')*pinvC*U';
%mw = H*v(:);

% All the cool kids are doing it this equivalent way;
diagC = diag(C,3);
diagCi = 1./diagC;
mw = X'\[zeros([3 1]); diagCi.*(U'*v(:))];

mom = reshape(mw(1:M*N), [N M])*diag(ht);
w = mw(M*N+1:end);

mwf = @(vel) subsref(X'\[zeros([3 1]); diagCi.*(U'*vel(:))], ...
                     struct('type', '()', 'subs', {{1:(M*N)}}));
v2p = @(vel) reshape(mwf(vel), [N M])*diag(ht);

pscl = @(p) subsref(p*diag(1/ht), struct('type', '()', 'subs', {{':'}}));
p2v = @(p) U*(diagC.*subsref(X'*[pscl(p); zeros([3 1])], ...
                      struct('type', '()', 'subs', {{4:M*N+3}})));

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
