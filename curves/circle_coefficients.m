function[k,s] = circle_coefficients(z)
% circle_coefficients -- Determines curvature coefficients for a circle
%
% [k,s] = circle_coefficients(z)
%
%     Given a 3-vector with complex entries z, this function determines the
%     2-vectors k and s such that p(t) = k(1) + k(2)*t is a primitive of the
%     curvature k(2) and satisfies:
%
%     z(2) = z(1) + \int_0^s(1) exp(i*p(t)) dt
%     z(3) = z(1) + \int_0^s(2) exp(i*p(t)) dt

z = z - z(1);

x = real(z);
y = imag(z);
a = abs(z);

k(1) = (y(2)*a(3)^2 - y(3)*a(2)^2)/(x(2)*a(3)^2 - x(3)*a(2)^2);
k(1) = atan2(y(2)*a(3)^2 - y(3)*a(2)^2, x(2)*a(3)^2 - x(3)*a(2)^2);

k(2) = -2*(x(2)*sin(k(1)) - y(2)*cos(k(1)))/a(2)^2;

tol = 1e-13;

if abs(k(2))<1e-13;  % Then it's still a line
  s(1) = abs(z(2)-z(1));
  s(2) = abs(z(3) - z(1));
else
  s(1) = atan2(k(2)*x(2)+sin(k(1)), -k(2)*y(2) + cos(k(1)));
  s(1) = (s(1) - k(1))/k(2);

  s(2) = atan2(k(2)*x(3)+sin(k(1)), -k(2)*y(3) + cos(k(1)));
  s(2) = (s(2) - k(1))/k(2);
end

s = mod(s, abs(2*pi/k(2)));  % s should be monotone

k = k(:);
s = s(:);
