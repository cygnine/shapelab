function[w, z, self] = calculcate_moebius_alignment(self, w, z)
% calculate_moebius_alignment -- Calculate and apply final Moebius map alignments for conformal maps
%
% [w, z, self] = moebius_alignment(self, w, z)
%
%     These final Moebius maps do four things:
%       1.) Connect user-specified shape-exterior points. (default: maps
%           infinity outside shape to infinity outside disc) This is a self-map
%           on H.
%           (self.moebius_maps.exterior_terminal)
%       2.) Connect user-specified shape-interior points. (default: do nothing)
%           This is a self-map on H.
%           (self.moebius_maps.interior_terminal)
%       3.) Rotates both the interior and exterior welding map points so that
%           the first sample point has angle 0. These are self-maps on D.
%           (Already built into maps 1 and 2, but separately accessible in
%            self.moebius_maps.interior_rotation and 
%            self.moebius_maps.exterior_rotation)
%       4.) Maps the result of all these points to the disc. 

%persistent moebius maps wrap
persistent wrap
if isempty(wrap)
  from labtools import interval_wrap as wrap
end

% Always specify the exterior map
map_to = 1/self.inf_image;

z_out_on_disc = self.moebius_maps.invert_D(...
                self.moebius_maps.H_to_D(z(end-2)));
%z_out_on_disc = moebius(z_n(end-2), maps.invert_D*maps.H_to_D);

maps.map_z_out_on_disc_to_0 = ...
     MoebiusMap([abs(z_out_on_disc)/z_out_on_disc, -abs(z_out_on_disc); ...
                 -conj(z_out_on_disc),             1]);

if abs(map_to)<(10*eps)
  maps.map_0_on_disc_to_goal = MoebiusMap(eye(2));
else
  maps.map_0_on_disc_to_goal = ...
     MoebiusMap([1,             abs(map_to); ...
                conj(map_to),  abs(map_to)/map_to]);
end

% To see what it does, follow it from the end
bigmap = maps.map_0_on_disc_to_goal.compose(...
             maps.map_z_out_on_disc_to_0.compose(...
                 self.moebius_maps.invert_D.compose(...
                     self.moebius_maps.H_to_D)));

% moebius(z_n(1), bigmap) is the location of the first sample point on the
% boundary of the disc. We want to rotate it to have angle 0.
exterior_rotation = -angle(bigmap(z(1)));

self.moebius_maps.exterior_rotation = MoebiusMap([exp(i*exterior_rotation), 0; ...
                                                 0,                        1]);

% Again, follow it from the end
self.moebius_maps.exterior_terminal = self.moebius_maps.D_to_H.compose(...
                                       self.moebius_maps.invert_D.compose(...
                                        self.moebius_maps.exterior_rotation.compose(...
                                         bigmap)));

% The above maps the real line to the real line...therefore normalize it to get
% rid of imaginary (machine eps crap) stuff
% God what I wouldn't give for a *= operator
self.moebius_maps.exterior_terminal = MoebiusMap(real(self.moebius_maps.exterior_terminal.H/...
                                       max(self.moebius_maps.exterior_terminal.H(:))));
  
z_in_on_disc = self.moebius_maps.H_to_D(z(end-1));

if isempty(self.center_location) & isempty(self.center_image);
  maps.map_z_in_on_disc_to_0 = MoebiusMap(eye(2));
else
  % TODO:
  % This assumes that center_image = 0.
  maps.map_z_in_on_disc_to_0 = ...
     MoebiusMap([abs(z_in_on_disc)/z_in_on_disc, -abs(z_in_on_disc); ...
                -conj(z_in_on_disc),            1]);
end

bigmap = maps.map_z_in_on_disc_to_0.compose(self.moebius_maps.H_to_D);

interior_rotation = -angle(bigmap(z(1)));

self.moebius_maps.interior_rotation = ...
       MoebiusMap([exp(i*interior_rotation), 0; ...
                   0,                        1]);

self.moebius_maps.interior_terminal = self.moebius_maps.D_to_H.compose(...
                                        self.moebius_maps.interior_rotation.compose(...
                                          bigmap));

self.moebius_maps.interior_terminal = MoebiusMap(...
      real(self.moebius_maps.interior_terminal.H/...
           max(self.moebius_maps.interior_terminal.H(:))));

zinf = z(end-2);
w = self.moebius_maps.interior_terminal(w);
z = self.moebius_maps.exterior_terminal(z);

% Update derivative
%H = self.moebius_maps.exterior_terminal.H;
%self.derivative_at_inf = self.derivative_at_inf*...
%  det(H)/(H(2,1)*zinf + H(2,2))^2;
self.derivative_at_inf = self.derivative_at_inf.*...
  self.moebius_maps.exterior_terminal.derivative(zinf);

% Before mapping to unit circle, also save fingerprint values
%self.interior_vertices = wrap(atan2(2*w(1:end-3), 1-w(1:end-3).^2), [0, 2*pi]);
%self.exterior_vertices = wrap(atan2(2*z(1:end-3), 1-z(1:end-3).^2), [0, 2*pi]);

wH = w;

% Finally, map to unit circle
w = self.moebius_maps.H_to_D(w);
z = self.moebius_maps.H_to_D(z);
% Compute final derivative: inf should now be at -i.
self.derivative_at_inf = -2*i*self.derivative_at_inf;

% Ok, remember we're doing M \circ \Psi \circ M^{-1}. Now do final M map taking
% Inf to -i, evaluate derivative there. Happy day: we're already at -i.
%self.derivative_at_inf = self.derivative_at_inf;

% Now rotate vertices so that the derivative is +real.
ang = angle(self.derivative_at_inf);
rotation = exp(i*(pi + ang));
%disp(pi + ang)

% finite difference approximation
ang2 = angle(diff(z(end-4:end-3)));
disp(ang2)
rotation = exp(i*-ang2);
%disp(ang2)

% Final rotation to fix derivative as real-valued
z = z*rotation;
self.moebius_maps.exterior_rotation = MoebiusMap([rotation 0; 0 1]);

% No exterior rotation
%self.moebius_maps.exterior_rotation = MoebiusMap(eye(2));

% Must reset interior rotation:
self.moebius_maps.interior_rotation = MoebiusMap(eye(2));
% Make real as usual
%self.moebius_maps.exterior_terminal = MoebiusMap(...
%    real(self.moebius_maps.exterior_terminal.H./...
%       max(self.moebius_maps.exterior_terminal.H(:))));

zH = real(self.moebius_maps.D_to_H(z));

self.interior_vertices = wrap(atan2(2*wH(1:self.N), 1-wH(1:self.N).^2), [0, 2*pi]);
self.exterior_vertices = wrap(atan2(2*zH(1:self.N), 1-zH(1:self.N).^2), [0, 2*pi]);

% Fix machine eps crap:
self.interior_vertices(1) = 0;
w(1) = 1;

self.interior_disc_vertices = w(1:self.N);
self.exterior_disc_vertices = z(1:self.N);

% Now unfortunately we want to make 0 match up with 0, so we have to figure out
% where on the interior this ends up.
temp = [6.2 linspace(0, 0.1, 2)];
ext0_on_int = self.interpolate([], temp);
% Why the F*** do I have to do the above with more than one point!?!?

% Rotate interior with this much
rotation = exp(-i*ext0_on_int(2));

self.moebius_maps.interior_rotation = MoebiusMap([rotation 0; 0 1]);

%self.moebius_maps.interior_terminal = self.moebius_maps.D_to_H.compose(...
%                  MoebiusMap([rotation 0; 0 1]).compose(...
%                     self.moebius_maps.H_to_D.compose(...
%                       self.moebius_maps.interior_terminal)));

w = rotation*w;
wH = real(self.moebius_maps.D_to_H(w));

self.interior_disc_vertices = w(1:self.N);
self.interior_vertices = wrap(atan2(2*wH(1:self.N), 1-wH(1:self.N).^2), [0, 2*pi]);
