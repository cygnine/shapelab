function[mu, nu] = conjugate_points(z1, z2)
% conjugate_points -- Finds conjugate points on  S^2
%
% [mu, nu] = conjugate_points(z1, z2)
%
%     Assuming that both z1 and z2 are different from 0, this function
%     identifies the geodesic on S2 passing through z1, z2, and 0, and returns
%     the two points mu and nu and are unchanged by assuming z1 and z2 to be
%     opposite ends of S^2.
%
%     This operation is unique -- obviously the choice of 0 as the third
%     point is arbitrary.

persistent moebius_inverse
if isempty(moebius_inverse)
  from shapelab.common import moebius_inverse
end

if isinf(z1) | isinf(z2)
  if isfinite(z1)  % Elegent boolean indexing doesn't do 0*Inf. Boo.
    A = z1;
  else
    A = z2;
  end
  %A = isfinite(z1)*z1 + isfinite(z2)*z2; % pick whichever is finite

  c = 2*A/(1-i);
  b = i*c;
  H = [1, b; -1, c];

  mu = moebius_inverse(0, H);
  nu = moebius_inverse(Inf, H);
else
  error('Oops, this hasn''t been coded since I haven''t needed it.');
end
