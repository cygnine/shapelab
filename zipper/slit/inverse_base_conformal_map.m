function[w] = inverse_base_conformal_map(z,a,varargin)
% [v,w] = inverse_base_conformal_map(z,a,{point_id=zeros(size(z)),cut_magnitude=abs(a))
%
%     Referring to [1], evaluates the inverse of the conformal mapping function
%     f_a, which is a basic building block for the `slit algorithm'
%     conformal mapping technique. The complex, scalar parameter is a, and the
%     (array) complex input is z. 
%
%     The optional input point_id has the same size as z and each entry takes
%     three possible values:
%     0: The point is some point in \mathbb{H}\backslash\mathbb{R}.
%     1: The point is located on \mathbb{R}
%     2: The point is located on \mathbb{R}, inside [p-1,p]. (See code for
%        definition of p.)
%
%  [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%  mapping", 2006.

persistent input_schema zipup
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import oblique_slit_zipup as zipup
end

opt = input_schema({'point_id'}, ...
      {zeros(size(z),'int8')}, [],varargin{:});

assert(length(a)==1, 'Error: not coded for vector-valued parameter a');

w = zipup(z,a,'point_id',opt.point_id);
