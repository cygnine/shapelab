function[dphi] = derivative(self,theta)
% derivative -- Derivative evaluation for Zipper-type welds
%
% dphi = derivative(self, theta)
%
%     Evaluates the derivative of the weld phi(theta) at the locations theta.
%     phi is the mapping from the 'interior' of [0,2*pi) to the 'exterior' of
%     [0, 2*pi).

% Since the weld has a (strictly) positive derivative, one can show that all we
% need is to evaluate the modulus of the derivative of Phi, where Phi is the map
% from S^1 to S^1. All other factors are angular and can be ignored. In
% practice, we'll use the real part of the logs to make this happen.

z = exp(i*theta);

tol = 1e-6; % maximum 'resolvable' distance between vertices and points.

% "Log DPhi"
Ldphi = zeros(size(z));

% self.inverse_moebius_alignment:
Ldphi = Ldphi + real(log(self.moebius_maps.interior_rotation.inv.derivative(z)));
z = self.moebius_maps.interior_rotation.inv(z);

Ldphi = Ldphi + real(log(self.moebius_maps.D_to_H.inv.derivative(z)));
z = self.moebius_maps.D_to_H(z);

% Fix real-valued crap:
z = real(z);
Ldphi = Ldphi + real(log(self.moebius_maps.interior_terminal.inv.derivative(z)));
z = self.moebius_maps.interior_terminal.inv(z);

% self.inverse_terminal_map: points not at 0
nzflags = abs(z) < tol;
zflags = not(nzflags);

sgn = -sign(self.winding_number);
z(zflags) = real(csqrt(sgn*z(zflags), 1/2, 'cut_bias', sgn>0));
Ldphi(zflags) = Ldphi(zflags) + real(log(sgn/2./z(zflags)));

Ldphi(zflags) = Ldphi(zflags) + self.moebius_maps.terminal_map.inv.derivative(z(zflags));
z(zflags) = self.moebius_maps.terminal_map.inv(z(zflags));

% self.inverse_terminal_map: points near 0: deal with N_teeth as well.
z(nzflags) = i*self.tooth_length;
Ldphi(nzflags) = Ldphi(nzflags) + log(1/(2*self.tooth_length));
Ldphi(nzflags) = Ldphi(nzflags) + real(log(self.moebius_maps.tooth_maps{end}.inv.derivative(i*self.tooth_length)));

% Still need to do N_teeth map for z(zflags)
