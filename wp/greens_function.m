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
%
%         G''' = -sin(theta)*log(2*(1-cos(theta)))

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
theta = wrap(theta, interval);
temp = (theta==0);
% These are flags for thetas that are close to 0 (mod 2 pi):
%temp = (abs(temp) < funeps) | (abs(temp - 2*pi) < funeps);

switch d
case 0
  stheta = sin(theta/2);
  G = 2*log(stheta) + log(4) - 3/2;
  G(temp) = 0;
  G = 1/2 + 2*stheta.^2.*G;

  %G = (1-cos(theta)).*log(2*(1-cos(theta))) + 3/2*cos(theta) - 1;
  %G(temp) = 1/2;
case 1
  G = sin(theta).*(2*log(sin(theta/2)) + log(4) - 1/2);

  %G = sin(theta).*log(2*(1-cos(theta))) - 1/2*sin(theta);
  G(temp) = 0;
case 2
  G = 2*cos(theta).*(log(sin(theta/2)) + log(2) + 1/4) + 1;

  %G = cos(theta).*log(2*(1-cos(theta))) + 1/2*cos(theta) + 1;
  G(temp) = 0;
case 3
  G = cos(theta/2).*cos(theta)./sin(theta/2) - sin(theta).*(1/2 + log(4) + 2*log(sin(theta/2)));

  %G = cot(theta).*(1 + cos(theta)) - sin(theta).*(1/2 + log(2*(1 - cos(theta))));
  G(temp) = 0;
case 4
  G = -1./(2*sin(theta/2).^2) - cos(theta).*(2*log(sin(theta/2)) + log(4) + 5/2) - 1;
  G(temp) = 0;
case 5
  stheta = sin(theta/2);
  G = sin(theta).*(2*log(stheta) + log(4) + 5/2 + 1./(2.*stheta.^2).*(1./(2*stheta.^2) - cos(theta)));
case 6
  stheta = sin(theta/2);
  G = cos(theta).*(2*log(stheta) + log(4) + 9/2 + (cos(theta) - 4)./(4*stheta.^4));
  G(temp) = 0;
otherwise
  error('Only implemented for d = 0,1,2,3,4,5,6');
end
