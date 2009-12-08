function[z] = polar_interpolate(z1, z2, N)
% polar_interpolate -- interpolates two points using a polar arc
%
% z = polar_interpolate(z1, z2, N)
% 
%     Interpolates the points z1 and z2 using a 'linspace' in polar coordinates,
%     meaning that the angles and radii of the output z are linearly spaced. It
%     does this in the 'shortest' way possible in the sense that it picks the
%     angular path with angular change no greater than pi.

r1 = abs(z1);
r2 = abs(z2);
a1 = mod(angle(z1), 2*pi);
a2 = mod(angle(z2), 2*pi);

if abs(a1-a2)<=pi
  z = linspace(r1, r2, N).*exp(i*linspace(a1, a2, N));
else
  if a2-a1>0
    a2 = a2 - 2*pi;
  else
    a2 = a2 + 2*pi;
  end
  z = linspace(r1, r2, N).*exp(i*linspace(a1, a2, N));
end
