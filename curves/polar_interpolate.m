function[z] = polar_interpolate(z1, z2, N)
% polar_interpolate -- interpolates two points using a polar arc
%
% z = polar_interpolate(z1, z2, N)
% 
%     Interpolates the points z1 and z2 using a 'linspace' in polar coordinates,
%     meaning that the angles and radii of the output z are linearly spaced.

r1 = abs(z1);
r2 = abs(z2);
a1 = mod(angle(z1), 2*pi);
a2 = mod(angle(z2), 2*pi);

z = linspace(r1, r2, N).*exp(i*linspace(a1, a2, N));
