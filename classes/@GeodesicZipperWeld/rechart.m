function[self] = rechart(self, varargin)
% rechart -- Recharts the welding map via a Moebius map applied to the interior
%
% [self] = rechart(self, interior_theta, exterior_theta)
%
%     If interior_theta and exterior_theta are length-3 vectors, then this
%     recharts the fingerprint so that the welding map maps interior_theta to
%     exterior_theta.
%
% [self] = rechart(self, center_location)
%
%     Effectively changes the protected property 'center_location' to the given
%     (scalar, complex-valued) location.

switch length(varargin)
case 1
  error('Not yet implemented');
case 2

  if not((length(varargin{1})==3) && (length(varargin{2})==3))
    error('Both inputs must be length-three real-valued vectors');
  end

  % Find where points are, and where we want them to be
  z_current = exp(i*self.interpolate([], varargin{2}));
  z_new = exp(i*varargin{1});

  M = MoebiusMap(z_new, z_current);
  M2 = self.moebius_maps.D_to_H.compose(...
            M.compose(...
                self.moebius_maps.H_to_D));
  M2 = MoebiusMap(inv(real(M2.H)));

  self.moebius_maps.interior_terminal = M2.compose(self.moebius_maps.interior_terminal);

  Mi = inv(M);
  % Now apply map to data in instance
  self.interior_disc_vertices = Mi(self.interior_disc_vertices);
  self.interior_vertices = angle(self.interior_disc_vertices);

  if abs(self.interior_vertices(1) - self.interior_interval(1)) < self.theta_tol
    self.interior_vertices(1) = self.interior_interval(1);
  end

  self.interior_interval = self.interior_vertices(1);

otherwise
  error('Method GeodesicZipperWeld.rechart accepts either one or two inputs');
end
