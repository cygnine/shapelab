function[c] = circle_grid(varargin)
% circle_grid -- Generates a grid of non-intersecting circles on the unit disc
%
% c = circle_grid({r=0.1,n=20})
%
%     This function generates a grid of tangential but non-intersecting circles
%     on the unit disc. The *approximate* radius is r. This is NOT circle packing; no
%     optimization is performed. It's a stupid radial stepping algorithm.
%     Simple, fast, and completely non-optimal. The matrix c that is returned
%     has size n x M, where n is the number of points on each circle and M is
%     the number of circles (determined at runtime). The samples are generated
%     as complex numbers.
%
%     Example:
%
%        c = circle_grid('r', 0.05, 'n', 35);
%        plot(c, 'k'); axis equal;

persistent strict_inputs bisection
if isempty(strict_inputs)
  from labtools import strict_inputs
  from labtools.rootfind import bisection
end

opt = strict_inputs({'r', 'n'}, {0.1, 20}, [], varargin{:});
r = opt.r;

M = 0;  % The total number of generated circles so far
R = 1;  % The current radius

% Here is a representative circle centered at 0 with radius 1:
circ = exp(i*linspace(0, 2*pi, opt.n).');

c = zeros([opt.n 0]);

while R>0
  if R<=opt.r % Then just add one circle and terminate
    c = [c R*circ];
    break
  end
  % Find the number of circles at the current radius

  thetafun = @(s) 2*atan(s/sqrt(R^2 - 2*R*s));

  dtheta = thetafun(r); % This is the angular span of a single circle of radius r
  M_current = ceil(2*pi/(dtheta)); % round up
  dtheta = 2*pi/(M_current);

  % Now find the r that creates angular span dtheta:
  r_current = bisection(0, r, thetafun, 'F', dtheta);

  % Now create the circles
  Mc = size(c, 2);
  c = [c zeros([opt.n M_current])];

  % Meh, could vectorize this...too lazy, not really worth it.
  for q = 1:M_current
    circ_center = (R-r_current)*exp(i*(q-1)*dtheta);
    c(:,Mc+q) = circ_center + r_current*circ;
  end

  R = R - 2*r_current;
end
