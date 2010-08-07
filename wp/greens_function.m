function[G] = greens_function(theta)
% greens_function -- Evaluates the Green's function of operator L
%
% G = greens_function(theta)
%
%     Evaluates the Green's function at the locations theta:
%
%         G = (1 - cos(theta)) * log(2*(1-cos(theta))) + 3/2*cos(theta) - 1

persistent wrap interval
if isempty(wrap)
  from labtools import interval_wrap as wrap
  interval = [0, 2*pi];
end

G = (1-cos(theta)).*log(2*(1-cos(theta))) + 3/2*cos(theta) - 1;
funeps = 1e-12;

temp = wrap(theta, interval);
G(abs(temp)<funeps) = 1/2;
G(abs(temp-2*pi)<funeps) = 1/2;
