function[G] = greens_function_derivative(theta)
% greens_function_derivative -- Evaluates derivative of the Green's function of operator L
%
% G = greens_function_derivative(theta)
%
%     Evaluates the Green's function derivative at the locations theta:
%
%         G = (1 - cos(theta)) * log(2*(1-cos(theta))) + 3/2*cos(theta) - 1
%
%         G'= sin(theta)*log(2*(1-cos(theta))) - 1/2*sin(theta)

persistent wrap interval
if isempty(wrap)
  from labtools import interval_wrap as wrap
  interval = [0, 2*pi];
end

G = sin(theta).*log(2*(1-cos(theta))) - 1/2*sin(theta);
funeps = 1e-12;

temp = wrap(theta, interval);
G(abs(temp)<funeps) = 0;
G(abs(temp-2*pi)<funeps) = 0;
