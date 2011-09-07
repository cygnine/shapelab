function[w, z, self] = calculate_terminal_map(self, w, z)
% calculate_terminal_map -- The terminal map for zipper-type algorithms
%
% [w, z, self] = calculate_terminal_map(self, w, z)
%
%     Determines and performs terminal maps for zipper-like algorithms.

%persistent moebius
%if isempty(moebius)
%  from shapelab.common import moebius
%end

self.a_array(end) = 1/z(1);

self.N_teeth = length(self.a_array)-1;

self.moebius_maps.terminal_map = MoebiusMap([1                      0; ...
                                             -self.a_array(end)  1]);

zinf = z(end-2);
w = self.moebius_maps.terminal_map(w);
z = self.moebius_maps.terminal_map(z);

%H = self.moebius_maps.terminal_map.H;
%self.derivative_at_inf = self.derivative_at_inf*det(H)/(H(2,1)*zinf + H(2,2))^2;
self.derivative_at_inf = self.derivative_at_inf*self.moebius_maps.terminal_map.derivative(zinf);

% Update derivative at inf
%self.derivative_at_inf = self.derivative_at_inf*...
%    (-2*sign(self.winding_number)*self.moebius_maps.terminal_map(zinf))/...
%    (-self.a_array(end)*zinf + 1)^2;


zinf = z(end-2);
w = -sign(self.winding_number)*w.^2;
z = -sign(self.winding_number)*z.^2;

zinf = self.moebius_maps.terminal_map(zinf);
self.derivative_at_inf = self.derivative_at_inf*...
   -2*sign(self.winding_number)*zinf;
