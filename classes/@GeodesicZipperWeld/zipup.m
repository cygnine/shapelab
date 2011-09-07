function[z_int, z_ext, self] = zipup(self, tooth, z_int, z_ext)
% zipup -- Zips up one tooth using the geodesic algorithm
%
% [z_int, z_ext, self] = zipup(self, tooth, z_int, z_ext);

persistent zipup
if isempty(zipup)
  from shapelab.loewner.solutions import normal_linear_slit_zip as zipup
end

zinf = z_ext(end);

% Things should already be symmetric, so do the zipup
z_int = zipup(i*self.tooth_length,z_int); 
z_ext = zipup(i*self.tooth_length,z_ext); 

self.derivative_at_inf = self.derivative_at_inf*zinf/z_ext(end);

if tooth > 1
  sgn = -sign(self.winding_number);
  % Now compute the next Moebius Map that takes z_int(tooth) and z_ext(tooth)
  % to symmetric locations.
  z1 = [sgn*self.tooth_length, 0, -sgn*self.tooth_length];
  z2 = [z_int(tooth), 0, z_ext(tooth)];

else
  % We just want to map z_int(1)==z_ext(1) to Inf
  % Also need to make sure self.inf_image winds up at 1
  z1 = [0, i, Inf];
  z2 = [0, z_ext(end), z_int(1)];
end

self.moebius_maps.tooth_maps{tooth} = MoebiusMap(z1, z2);

zinf = z_ext(end);
z_int = self.moebius_maps.tooth_maps{tooth}.inv(z_int);
z_ext = self.moebius_maps.tooth_maps{tooth}.inv(z_ext);

%H = self.moebius_maps.tooth_maps{tooth}.inv.H;
%self.derivative_at_inf = self.derivative_at_inf*...
%   det(H)/(H(2,1)*zinf + H(2,2))^2;
self.derivative_at_inf = self.derivative_at_inf*...
    self.moebius_maps.tooth_maps{tooth}.inv.derivative(zinf);
