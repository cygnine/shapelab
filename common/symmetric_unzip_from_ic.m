function[v,w] = symmetric_unzip_from_ic(z,c,varargin)
% [v,w] = symmetric_unzip_from_ic(z,c,{point_id=zeros(size(z))})
%
%     Implements the rather unorthodox mapping sqrt(z^2 + c^2), for some
%     real-valued c. This mapping takes the line segment (0,0) -- (0,c) and
%     'unzips' it into two line segments (-c,0) -- (0,0) and (0,0) -- (c,0).
%     Nothing is unusual about this map except on the x-axis: the sign of the
%     input is preserved. 
%
%     The optional input point_id has the same size as z and each entry takes
%     three possible values:
%     0: The point is some point in \mathbb{H}\backslash\{\gamma}, where \gamma
%        is the line segment (0,0) -- (0,c).
%     1: The point is located on \gamma
%     2: The point is located on \partial\mathbb{H} = \mathbb{R}
%
%     The behavior of this map is as follows for each of these possibilities:
%     0: The map takes \mathbb{H}\backslash\{\gamma} \rightarrow
%        \mathbb{H}\backslash\{\gamma}, i.e. no fancy stuff is required.
%     1: The map 'unzips' \gamma into two segments: (-c,0) -- (0,0) and (0,0) --
%        (c,0). Fort this reason, there are two outputs: v,w. v corresponds to
%        the (-c,0) segment and w corresponds to the (c,0) segment. 
%     2: The map is standard, but preserves the sign of z. 
%
%     The operations 1 and 2 are 'one-dimensional'; the operation on 0 is a
%     complex-number operation.

global handles;
%opt = handles.common.InputSchema({'cut_bias','boundary'}, {false,false}, [], varargin{:});
opt = handles.common.InputSchema({'point_id'}, {zeros(size(z),'int8')}, [], varargin{:});
csqrt = handles.shapelab.common.positive_angle_square_root;

interior = opt.point_id==0;
gamma = opt.point_id==1;
boundary = opt.point_id==2;

v = zeros(size(z));
w = zeros(size(z));

% These points are supposed to be on the real line
z(boundary) = real(z(boundary));

% These points are supposed to be purely imaginary
z(gamma) = i*imag(z(gamma));

v(interior) = csqrt(z(interior).^2 + c^2);
v(boundary) = sign(z(boundary)).*sqrt(z(boundary).^2 + c^2);
w = v;

v(gamma) = -sqrt(c^2 - imag(z(gamma)).^2);
w(gamma) = -v(gamma);
