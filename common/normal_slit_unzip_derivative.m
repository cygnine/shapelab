function[v,w] = normal_slit_unzip_derivative(z,c,varargin)
% [v,w] = normal_slit_unzip_derivative(z,c,{point_id=zeros(size(z))})
%
%     Implements the derivative of normal_slit_unzip.
%
%     The optional input point_id has the same size as z and each entry takes
%     three possible values:
%     0: The point is some point in \mathbb{H}\backslash\{\gamma}, where \gamma
%        is the line segment (0,0) -- (0,c).
%     1: The point is located on \gamma
%     2: The point is located on \partial\mathbb{H} = \mathbb{R}
%
%     TODO: the derivative on gamma (point_id=1) is super-sketchy.

persistent input_schema csqrt
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import positive_angle_square_root as csqrt
end

opt = input_schema({'point_id'}, {zeros(size(z),'int8')}, [], varargin{:});

interior = opt.point_id==0;
gamma = opt.point_id==1;
boundary = opt.point_id==2;

v = zeros(size(z));
w = zeros(size(z));

% These points are supposed to be on the real line
z(boundary) = real(z(boundary));

% These points are supposed to be purely imaginary
z(gamma) = i*imag(z(gamma));

%v(interior) = csqrt(z(interior).^2 + c^2);
v(interior) = z(interior)./(csqrt(z(interior).^2 + c^2));

%v(boundary) = sign(z(boundary)).*sqrt(z(boundary).^2 + c^2);
v(boundary) = sign(z(boundary)).*z(boundary)./sqrt(z(boundary).^2 + c^2);

w = v;

%v(gamma) = -sqrt(c^2 - imag(z(gamma)).^2);
v(gamma) = -z(gamma)./sqrt(c^2 - imag(z(gamma)).^2);
%w(gamma) = -v(gamma);
w(gamma) = -v(gamma);
