function[z] = symmetric_zipup_to_ic(w,c,varargin)
% [z] = symmetric_zipup_to_ic(w,c,point_id=zeros(size(z)))
%
%     Implements the inverse of symmetric_unzip_from_ic. I.e., the function sqrt(z^2 -
%     c^2) for a real-valued c. As in symmetric_unzip_from_ic, the optional input point_id
%     determines how the map behaves. The possible values of point_id are:
%
%     0: The point is somewhere in \mathbb{H}\backslash{\mathbb{R}}. Action of
%        the map is normal.
%     1: The point is on \mathbb{R}. 
%     2: The point is on the \mathbb{R} between [-c,c]. Although this is an
%     easily-testable case of 1, you may want to force this behavior in the case
%     of machine-epsilon crap.
%
%     The number of possibilities are reduced compared to symmetric_unzip_from_ic because
%     there are really only two regions where care must be taken.

persistent csqrt input_schema
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import positive_angle_square_root as csqrt
end

opt = input_schema({'point_id'}, {zeros(size(w))}, [], varargin{:});

interior = opt.point_id==0;
rline = opt.point_id==1;
gamma = opt.point_id==2;

% The 'interior' case
z(interior) = csqrt(w(interior).^2 - c^2);

% The problem cases
z1 = real(w(rline));
z2 = real(w(gamma));
% Let's do some massaging. 'gamma' means on the pre-image of (0,0) -- (0,c)
z1gamma = abs(z1)<=c; % cases where gamma behavior is not explicit
z1(z1gamma) = i*sqrt(c^2 - z1(z1gamma).^2);
z1(~z1gamma) = sign(z1(~z1gamma)).*sqrt(z1(~z1gamma).^2 - c^2);
z(rline) = z1;

% If things are on gamma, straightforward
z2 = i*sqrt(c^2 - z2.^2);
z(gamma) = z2;
