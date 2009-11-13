function[z] = normal_slit_zipup_derivative(w,c,varargin)
% [z] = normal_slit_zipup_derivative(w,c,point_id=zeros(size(z)))
%
%     Implements the derivative of normal_slit_zipup, i.e. of the function sqrt(z^2 -
%     c^2) for a real-valued c. As in normal_slit_unzip, the optional input point_id
%     determines how the map behaves. The possible values of point_id are:
%
%     0: The point is somewhere in \mathbb{H}\backslash{\mathbb{R}}. Action of
%        the map is normal.
%     1: The point is on \mathbb{R}. 
%     2: The point is on the \mathbb{R} between [-c,c]. Although this is an
%        easily-testable case of 1, you may want to force this behavior in the case
%        of machine-epsilon crap.

persistent csqrt input_schema
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import positive_angle_square_root as csqrt
end

opt = input_schema({'point_id'}, {zeros(size(w))}, [], varargin{:});

interior = opt.point_id==0;
% Lump all cases when w is on the real line together:
rline = (opt.point_id==1) | (opt.point_id==2);

% The 'interior' case
z(interior) = w(interior)./csqrt(w(interior).^2 - c^2);

% I think for the real line things are still straightforward:
z(rline) = w(rline)./sqrt(w(rline) - c^2);
