function[z_int, z_ext, self] = calculate_inverse_moebius_alignment(self, z_int, z_ext);
% calculate_inverse_moebius_alignment -- Computes inverse moebius alignment
%
% [z_int, z_ext, self] = calculate_inverse_moebius_alignment(self, z_int, z_ext)

% Need to find the map on H -- first we want z_int(end-1), z_ext(end-1) to end up at
% 0, and we want z_int(1) and z_ext(1) to end up at Inf.

inner_rotation = conj(z_int(1));
outer_rotation = conj(z_ext(1));

% First rotate appropriately to 0
%self.moebius_maps.interior_rotation = MoebiusMap([conj(inner_rotation) 0; 0 1]);
%self.moebius_maps.exterior_rotation = MoebiusMap([conj(outer_rotation) 0; 0 1]);
self.moebius_maps.interior_rotation = MoebiusMap(eye(2));
self.moebius_maps.exterior_rotation = MoebiusMap(eye(2));

z_int = self.moebius_maps.interior_rotation.inv(z_int);
z_ext = self.moebius_maps.exterior_rotation.inv(z_ext);

H = self.moebius_maps.exterior_rotation.inv.H;
self.derivative_at_inf = det(H);

z_int = self.moebius_maps.D_to_H(z_int);
z_ext = self.moebius_maps.D_to_H(z_ext);

z_int(1:self.N) = real(z_int(1:self.N));
z_ext(1:self.N) = real(z_ext(1:self.N));

%self.derivative_at_inf = -2*i*self.derivative_at_inf;
self.derivative_at_inf = i/2*self.derivative_at_inf;

% Ok, remember we're doing M \circ \Psi \circ M^{-1}. Now do final M map taking
% Inf to -i, evaluate derivative there. Happy day: we're already at -i.
%self.derivative_at_inf = self.derivative_at_inf;

% Now specify where things should go
z1 = [0, Inf, -self.tooth_length];
z2 = [z_int(self.N), z_int(1), z_int(self.N-1)];
temp = MoebiusMap(z1, z2);
self.moebius_maps.interior_terminal = MoebiusMap(temp.H/max(temp.H(:)));
self.moebius_maps.interior_terminal = MoebiusMap(...
  real(self.moebius_maps.interior_terminal.H));

z1 = [0, Inf, -self.tooth_length];
z2 = [z_ext(self.N), z_ext(1), z_ext(self.N-1)];
temp = MoebiusMap(z1, z2);
self.moebius_maps.exterior_terminal = MoebiusMap(temp.H/max(temp.H(:)));
self.moebius_maps.exterior_terminal = MoebiusMap(...
  real(self.moebius_maps.exterior_terminal.H));

zinf = z_ext(end);
z_int = self.moebius_maps.interior_terminal.inv(z_int);
z_ext = self.moebius_maps.exterior_terminal.inv(z_ext);

%H = self.moebius_maps.exterior_terminal.inv.H;
%self.derivative_at_inf = self.derivative_at_inf*...
%   det(H)/(H(2,1)*zinf + H(2,2))^2;
self.derivative_at_inf = self.derivative_at_inf*...
   self.moebius_maps.exterior_terminal.inv.derivative(zinf);
