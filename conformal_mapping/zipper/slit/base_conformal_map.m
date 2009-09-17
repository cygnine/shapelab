function[v,w] = base_conformal_map(z,a,varargin)
% [v,w] = base_conformal_map(z,a,{point_id=zeros(size(z)),cut_magnitude=abs(a)^2/imag(a))
%
%     Referring to [1], evaluates the conformal mapping function f_a, which is a
%     basic building block for the `slit' conformal mapping
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
%        is line segment (0,0) -- a.
%     1: The point is located on \gamma
%     2: The point is located on \partial\mathbb{H} = \mathbb{R}
%
%  [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%       mapping", 2006.

global handles;
opt = handles.common.input_schema({'point_id'}, ...
      {zeros(size(z),'int8')}, [],varargin{:});
%unzip = handles.shapelab.common.slit_unzip_from_a;
unzip = handles.shapelab.conformal_mapping.zipper.slit.slit_unzip_from_a;

assert(length(a)==1, 'Error: not coded for vector-valued parameter a');

[v,w] = unzip(z,a,'point_id', opt.point_id);
