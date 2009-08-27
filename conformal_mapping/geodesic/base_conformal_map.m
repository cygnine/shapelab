function[v,w] = base_conformal_map(z,a,varargin)
% [v,w] = base_conformal_map(z,a,{point_id=zeros(size(z)),cut_magnitude=abs(a)^2/imag(a))
%
%     Referring to [1], evaluates the conformal mapping function f_a, which is a
%     basic building block for the `geodesic algorithm' conformal mapping
%     technique. The complex, scalar parameter is a, and the (array) complex
%     input is z. 
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

global handles;
% Intermediate points: see [1]
b = abs(a)^2/real(a);
c = abs(a)^2/imag(a);
opt = handles.common.InputSchema({'point_id','cut_magnitude'}, ...
      {zeros(size(z),'int8'), c}, [],varargin{:});
moebius = handles.shapelab.common.moebius;
unzip_at_c = handles.shapelab.common.z_unzip_at_ic;

assert(length(a)==1, 'Error: not coded for vector-valued parameter a');

factor = opt.cut_magnitude/c;

w = moebius(z, [factor 0;-1/b 1]);
[v,w] = unzip_at_c(w,opt.cut_magnitude,'point_id', opt.point_id);
