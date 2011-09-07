function[z_int, z_ext, self] = calculate_inverse_terminal_map(self, z_int, z_ext)
% calculate_inverse_terminal_map -- The inverse terminal map for zipper-type algorithms
%
% [z_int, z_ext, self] = calculate_inverse_terminal_map(self, z_int, z_ext)
%
%     Determines and performs inverse terminal maps for zipper-like algorithms.

persistent csqrt
if isempty(csqrt)
  from shapelab.common import positive_angle_exponential as csqrt
end

sgn = -sign(self.winding_number);

zinf = z_ext(end);

% Do sqrt business
z_int = csqrt(sgn*z_int, 1/2, 'cut_bias', sgn>0);
z_ext = csqrt(sgn*z_ext, 1/2, 'cut_bias', sgn<0);

%self.derivative_at_inf = self.derivative_at_inf*1/2*sgn/z_ext(end);
self.derivative_at_inf = self.derivative_at_inf*sqrt(sgn)*1/2/sqrt(zinf);

z_int(1:self.N) = real(z_int(1:self.N));
z_ext(1:self.N) = real(z_ext(1:self.N));

% Now z_int(end-2) and z_ext(end-2) are unsymmetric around 0...find the Moebius
% map that makes them symmetric.
z1 = [sgn*self.tooth_length, 0, -sgn*self.tooth_length];
z2 = [z_int(self.N-1), 0, z_ext(self.N-1)];
self.moebius_maps.terminal_map = MoebiusMap(z1, z2);

zinf = z_ext(end);

% Apply the map
z_int = self.moebius_maps.terminal_map.inv(z_int);
z_ext = self.moebius_maps.terminal_map.inv(z_ext);

%H = self.moebius_maps.terminal_map.inv.H;
%self.derivative_at_inf = self.derivative_at_inf*...
%   det(H)/(H(2,1)*zinf + H(2,2))^2;
self.derivative_at_inf = self.derivative_at_inf*...
  self.moebius_maps.terminal_map.inv.derivative(zinf);
