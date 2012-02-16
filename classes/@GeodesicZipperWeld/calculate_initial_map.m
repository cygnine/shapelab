function[z,self] = calculate_initial_map(self, z)
% calculate_initial_map -- The initial map for the geodesic zipper-type algorithm.
%
% z = calculate_initial_map(self, z)
%
%     Performs the initial map for zipper-type algorithms. 

self.moebius_maps.initial_map = MoebiusMap([1 -z(2); ...
                                            1 -z(1)]);

z([1:2 end]) = [Inf; 0; Inf];

z(3:end-1) = i*sqrt(self.moebius_maps.initial_map(z(3:end-1)));

% At this point, infinity has been mapped to +i. Compute derivative, modulo
% 1/inf^2 factor up to this point.
self.derivative_at_inf = i/2*(self.z(2) - self.z(1));

% Let's just see if we can get *most* of the map derivatives computed, modulo
% the maps to/from infinity
%z(end-4) = i;
%z(end-3) = i + 1e-4;
self.derivative_at_inf = 1;

% Ok, we're actually going to find the derivative of the map M \circ \Psi \circ
% M^{-1}, where M is a Moebius map that takes Inf to -i (the D ---> H map). The
% reason is so that we can avoid dealing with infinity.

% Find map taking -i ---> Inf ---> +1 (emulates moebius_maps.initial_map)
%temp = self.moebius_maps.initial_map.compose(self.moebius_maps.H_to_D);
%H = temp.H;
%self.derivative_at_inf = det(temp.H)/(H(2,1)*(-i) + H(2,2))^2;
%self.derivative_at_inf = i/2*self.derivative_at_inf;
