function[z_int, z_ext, self] = calculate_inverse_initial_map(self, z_int, z_ext);
% calculate_inverse_initial_map -- The inverse initial map for zipper-type algorithms
%
% [z_int, z_ext, self] = calculate_inverse_initial_map(self, z_int, z_ext)
%
%     Performs the inverse of the initial map for zipper-type algorithms.

zinf = z_ext(end);

% First step is easy:
z_int = (-i*z_int).^2;
z_ext = (-i*z_ext).^2;

self.derivative_at_inf = self.derivative_at_inf*(-2*zinf);

% Now how to compute the initial Moebius Map? We have the images of one interior
% and one exterior point (center/inf). Make sure they map to
% self.center_location and self.inf_location using a map whose inverse is
% the two-parameter Moebius map
% 
%       z - z_0
% M =  ---------
%       z - z_1
%
% The numbers z_0 and z_1 are the original locations of z_int(1:2)==z_ext(1:2).
%
% By construction, z_ext(end) (the image of self.inf_image) should be at +1 now.
% No matter what, the map M will send it back to Inf.
% There is one degree of freedom:
z0 = 1; % Completely random

z1 = [self.inf_location, self.center_location, z0];
z2 = [z_ext(end), z_int(end), Inf];

self.moebius_maps.initial_map = MoebiusMap(z1, z2);

z_int = self.moebius_maps.initial_map.inv(z_int);
z_ext = self.moebius_maps.initial_map.inv(z_ext);

H = self.moebius_maps.initial_map.inv.H;
self.derivative_at_inf = self.derivative_at_inf*...
   det(H)/H(2,1)^2;

%disp(angle(det(H)/H(2,1)^2))
% Find map taking +1 ---> Inf ---> -i (emulates moebius_maps.terminal_map)
%temp = self.moebius_maps.D_to_H.compose(self.moebius_maps.terminal_map);
%H = temp.H;
%self.derivative_at_inf = det(temp.H)/(H(2,1)*(1) + H(2,2))^2;

% Now we need to rotate so that the derivative at inf is a positive real.
ang = angle(self.derivative_at_inf);
disp(ang);
%rotation = exp(i*ang);

% finite-difference approximation
ang2 = angle(diff(z_ext(self.N+1:self.N+2)));
%disp(ang2)
rotation = exp(-i*ang2);

temp = MoebiusMap([rotation 0; 0 1]);
z_int = temp(z_int);
z_ext = temp(z_ext);

%% Update initial map
temp = temp.compose(self.moebius_maps.initial_map.inv);
self.moebius_maps.initial_map = inv(temp);

% Now we have to shift the image of the center back to the center, and we again
% make sure |z(1)| = 1 -- completely random.
temp = MoebiusMap(1/abs(z_int(1))*[1 -z_int(end); 0 1]);
z_int = temp(z_int);
z_ext = temp(z_int);
temp = temp.compose(self.moebius_maps.initial_map.inv);
self.moebius_maps.initial_map = inv(temp);

self.z = 1/2*(z_int(1:self.N) + z_ext(1:self.N));
