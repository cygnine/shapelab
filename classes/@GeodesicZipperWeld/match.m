function[self] = match(self, other)
% match -- Performs a scale+shift map on the shape self
%
% self = match(self, other)
%
%     Given another GeodesicZipperWeld object "other", this function computes
%     the Moebius scale + shift that minimizes the least-squares difference
%     between points. It does this by treating the points on the "other" shape
%     as landmarks, and tries to match points on the current shape to them.

% We'll do N least-squares optimizations
N = length(self.z);
errors = zeros([N 1]);
maps = zeros([3 N]);

w = other.z;
other_theta_ext_offset = other.exterior_vertices - other.exterior_vertices(1);
M = length(w);
w = [real(w); imag(w)];
% We'll try to scale and shift self so that M points on self match with the
% above M points on other.

for q = 1:N;
  % Get M points on self starting with self.z(q)
  on_exterior_self = self.exterior_vertices(q) + other_theta_ext_offset;
  z = self.map_to_shape(exp(i*on_exterior_self), true([M 1]));

  % Find scale+shift that minimize least-squares here:
  A = [real(z) ones(size(z)) zeros(size(z)); ...
       imag(z) zeros(size(z)) ones(size(z))];
  maps(:,q) = A\w;

  errors(q) = norm(A*maps(:,q) - w);
end

[garbage, q] = min(errors);
% So we use scale and shift q
H = [maps(1,q) maps(2,q) + i*maps(3,q); ...
        0         1];

M = MoebiusMap(H);

self.z = M(self.z);
%self.moebius_maps.initial_map = M.inv.compose(self.moebius_maps.initial_map);
self.moebius_maps.initial_map = self.moebius_maps.initial_map.compose(M.inv);
