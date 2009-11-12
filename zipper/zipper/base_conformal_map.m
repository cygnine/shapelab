function[v,w] = base_conformal_map(z,c,a,varargin)
% [v,w] = base_conformal_map(z,c,a,{point_id=zeros(size(z)),cut_magnitude=abs(a)^2/imag(a))
%
%     Referring to [1], evaluates the conformal mapping function f_{c,a}, which is a
%     basic building block for the `zipper' conformal mapping
%     technique. The complex, scalar parameters are a and c, and the (array)
%     complex input is z. The parameter a is the 'tip' of the circular
%     arc, whereas c is between 0 and a on the arc.
%
%     The two outputs v and w refer to whether the points that lie on gamma (see
%     below) are unzipped to the negative x-axis or positive x-axis,
%     respectively.
%
%     The optional input point_id has the same size as z and each entry takes
%     three possible values:
%     0: The point is some point in \mathbb{H}\backslash\{\gamma}, where \gamma
%        is line segment (0,0) -- a.
%     1: The point is located on \gamma
%     2: The point is located on \partial\mathbb{H} = \mathbb{R}
%
%  [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%       mapping", 2006.

persistent input_schema unzip moebius
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.conformal_mapping.zipper.slit import slit_unzip_from_a as unzip
  from shapelab.common import moebius
end

opt = input_schema({'point_id'}, ...
      {zeros(size(z),'int8')}, [],varargin{:});

assert(length(a)==1, 'Error: not coded for vector-valued parameter a');
assert(length(c)==1, 'Error: not coded for vector-valued parameter c');

delta_real = real(a) - real(c);
if abs(delta_real)<1e-14
  d = a;  % Straight line segment...we *could* actually use geodesic here since
          % it's the same conformal map
  m = eye(2);
else
  b = (imag(c)*abs(a)^2 - imag(a)*abs(c)^2)/(real(a)*imag(c) - real(c)*imag(a));
  m = [1 0; -1/b 1];
  z = moebius(z,m);
  d = moebius(a,m);
end

[v,w] = unzip(z,d,'point_id', opt.point_id);
