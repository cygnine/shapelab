function[G] = greens_function(theta, d)
% greens_function -- Evaluates the Green's function of operator L
%
% G = greens_function(theta, [[d==0]])
%
%     Evaluates the d'th derivative of the Green's function at the locations theta:
%
%         G = (1 - cos(theta)) * log(2*(1-cos(theta))) + 3/2*cos(theta) - 1
%
%         G'= sin(theta)*log(2*(1-cos(theta))) - 1/2*sin(theta)
%         G'(0) = 0   (continuous limit)
%
%         G'' = cos(theta)*log(2*(1-cos(theta))) + 1/2*cos(theta) + 1, 
%         G''(0) = 0   (definition)

persistent wrap interval
if isempty(wrap)
  from labtools import interval_wrap as wrap
  interval = [0, 2*pi];
end

if nargin < 2
  d = 0;
end

% Flags for evaluating close to theta=0
funeps = 1.1e-8;
temp = wrap(theta, interval);
% These are flags for thetas that are close to 0 (mod 2 pi):
temp = (abs(temp) < funeps) | (abs(temp - 2*pi) < funeps);

switch d
case 0
  G = (1-cos(theta)).*log(2*(1-cos(theta))) + 3/2*cos(theta) - 1;
  G(temp) = 1/2;
  %G(abs(temp)<funeps) = 1/2;
  %G(abs(temp-2*pi)<funeps) = 1/2;
case 1
  G = sin(theta).*log(2*(1-cos(theta))) - 1/2*sin(theta);
  G(temp) = 0;
  %G(abs(temp)<funeps) = 0;
  %G(abs(temp-2*pi)<funeps) = 0;
case 2
  G = cos(theta).*log(2*(1-cos(theta))) + 1/2*cos(theta) + 1;
  G(temp) = 0;
  %G(abs(temp)<funeps) = 0;
  %G(abs(temp-2*pi)<funeps) = 0;
otherwise
  error('Only implemented for d = 0,1,2');
end
