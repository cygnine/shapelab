function[G] = greens_function_derivative(theta, d)
% greens_function_derivative -- Evaluates derivative of the Green's function of operator L
%
% G = greens_function_derivative(theta, [[d==1]])
%
%     Evaluates the d'th order Green's function derivative at the locations theta:
%
%         G = (1 - cos(theta)) * log(2*(1-cos(theta))) + 3/2*cos(theta) - 1
%
%         G'= sin(theta)*log(2*(1-cos(theta))) - 1/2*sin(theta)
%
%         G'' = cos(theta)*log(2*(1-cos(theta))) + 1/2*cos(theta) + 1, 
%         G''(0) = 0   (definition)

persistent wrap interval
if isempty(wrap)
  from labtools import interval_wrap as wrap
  interval = [0, 2*pi];
end

fprintf('This function is obsolete. Use shapelab.wp.greens_function instead\n');

if nargin < 2
  d = 1;
end

funeps = 1.1e-8;
temp = wrap(theta, interval);

switch d
case 1
  G = sin(theta).*log(2*(1-cos(theta))) - 1/2*sin(theta);

  G(abs(temp)<funeps) = 0;
  G(abs(temp-2*pi)<funeps) = 0;
case 2
  G = cos(theta).*log(2*(1-cos(theta))) + 1/2*cos(theta) + 1;
  G(abs(temp)<funeps) = 0;
  G(abs(temp-2*pi)<funeps) = 0;
otherwise
  error('Only implemented for d=1 or d=2');
end
