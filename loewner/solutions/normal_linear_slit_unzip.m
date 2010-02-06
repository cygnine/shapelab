function[w, w_slit_left, w_slit_right] = normal_linear_slit_unzip(a, z, varargin)
% normal_linear_slit_unzip -- Unzips a linear slit that is normal to the real axis
%
% [w, w_slit_left, w_slit_right] = normal_linear_slit_unzip(a, z, ...
%                                   {z_slit = []})
%
%     The exact solution to the Loewner equation for a linear slit that is
%     perpendicular to the axis at the location x = real(a).
%
%                   a
%                     x
%                     |  
%      z-plane        |                            w-plane
%                     |                   ----->     
%                     |
%                     |                                           f(a) = real(a)
%     ________________o________________           _________o-----x-----o________
%
%
%     The input vector z contains all points that are not on the slit; the map
%     in this case is simply implemented as there is no ambiguity. The optional
%     input z_slit contains points on the slit that are to be unzipped. No
%     testing is performed on z_slit to make sure it actually does lie on the
%     slit. If z_slit is given, then both the `left' and `right' images are
%     given in output.
%
%     If you just want to evaluate the map for points on the slit, set the
%     mandatory input z to the empty array ([]).

persistent cexp bsign
if isempty(cexp)
  from shapelab.common import positive_angle_exponential as cexp
  from labtools import biased_sign as bsign
end

if nargin>2
  z_slit = varargin{1};
  slit = true;
else
  slit = false;
end

w = zeros(size(z));
ra = real(a);  % The point on the real axis to unzip to 
ia = imag(a);

% All the regular points:
if not(isempty(z))
  reals = (imag(z)==0);
  nreals = not(reals);

  if any(reals)
    w(reals) = bsign(z(reals)-ra).*sqrt((z(reals)-ra).^2 + ia^2) + ra;
  end

  if any(nreals)
    w(nreals) = cexp((z(nreals) - ra).^2 + ia^2, 1/2) + ra;
  end
end

if slit
  w_slit_left = zeros(size(z_slit));
  w_slit_right = zeros(size(z_slit));

  discriminant = sqrt(-imag(z_slit).^2 + ia^2);
  w_slit_left = ra - discriminant;
  w_slit_right = ra + discriminant;
else
  w_slit_left = [];
  w_slit_right = [];
end
