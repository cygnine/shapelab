function[F] = curvature_function(k,s,z)
% curvature_function -- evaluates nonlinear function describing curvature
%
% F = curvature_function(k,s,z)
%
%   If length(z) == length(s) == length(k)+1, then this function evaluates returns
%   a vector of length 2*length(z). Let x = real(z) and y = imag(z). Then the
%   entries are:
%
%   F(1) = \int_0^s(1) cos(k(1) + k(2)*s + k(3)*s^2 + ...) ds - (x(2) - x(1))
%   F(2) = \int_0^s(1) sin(k(1) + k(2)*s + k(3)*s^2 + ...) ds - (y(2) - y(1))
%
%   F(3) = \int_s(1)^s(2) cos(k(1) + k(2)*s + k(3)*s^2 + ...) ds - (x(3) - x(2))
%   F(4) = \int_s(1)^s(2) sin(k(1) + k(2)*s + k(3)*s^2 + ...) ds - (y(3) - y(2))
%      .
%      .
%      .
%      .
%
%   Note: this function does not assume that z(1) = 0, but clearly from the
%   definition of the function evaluations, this is meant to be the case.

persistent trig_integral
if isempty(trig_integral)
  from shapelab.curves import trig_integral
end

N = length(s);
s = s(:);
F = zeros([2*N 1]);

cells = zeros([N 2]);
cells(:,1) = [0; s(1:N-1)];
cells(:,2) = s(:);

[cosint,sinint] = trig_integral(k,cells);

z = diff(z(:));
cosint = cosint - real(z);
sinint = sinint - imag(z);

temp = [cosint, sinint].';
F = temp(:);
