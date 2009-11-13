function[v,w] = base_map(z,a,varargin)
% [v,w] = base_map(z,a,{point_id=zeros(size(z)),cut_magnitude=abs(a)^2/imag(a))
%
%     Referring to [1], evaluates the conformal mapping function f_a, which is a
%     basic building block for the `geodesic algorithm' conformal mapping
%     technique. The complex, scalar parameter is a, and the (array) complex
%     input is z. 
%
%     The two outputs v and w refer to whether the points that lie on gamma (see
%     below) are unzipped to the negative x-axis or positive x-axis,
%     respectively.
%
%     The optional input point_id has the same size as z and each entry takes
%     three possible values:
%     0: The point is some point in \mathbb{H}\backslash\{\gamma}, where \gamma
%        is the circular arc (0,0) -- a, orthogonal to \mathbb{R}.
%     1: The point is located on \gamma
%     2: The point is located on \partial\mathbb{H} = \mathbb{R}
%
%     cut_magnitude refers to what position on the x-axis the point a ends up
%     at. The default is simply what the default is in Marshall's paper.
%
%  [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%       mapping", 2006.

persistent input_schema moebius unzip_at_c
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import moebius
  from shapelab.common import normal_slit_unzip as unzip_at_c
end

% Intermediate points: see [1]
b = abs(a)^2/real(a);
c = abs(a)^2/imag(a);
opt = input_schema({'point_id','cut_magnitude'}, ...
      {zeros(size(z),'int8'), c}, [],varargin{:});

assert(length(a)==1, 'Error: not coded for vector-valued parameter a');

factor = opt.cut_magnitude/c;

w = moebius(z, [factor 0;-1/b 1]);
[v,w] = unzip_at_c(w,opt.cut_magnitude,'point_id', opt.point_id);
