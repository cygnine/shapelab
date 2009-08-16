function[w] = z_unzip_at_ic(z,c,varargin)
% [w] = z_unzip_at_ic(z,c,{cut_bias=false,boundary=false})
%
%     Implements the rather unorthodox mapping sqrt(z^2 + c^2), for some
%     real-valued c. This mapping takes the line segment (0,0) -- (0,c) and
%     'unzips' it into two line segments (-c,0) -- (0,0) and (0,0) -- (c,0).
%     Nothing is unusual about this map except on the x-axis: the sign of the
%     input is preserved. 
%
%     The optional input cut_bias indicates which side of the x-axis the line
%     segment (0,0) -- (0,c) gets mapped to. true = (+) x-axis, false = (-)
%     x-axis.
%
%     If boundary is set to true, this function assumes that all inputs must lie
%     on the boundary of a shape and so there cannot be any inputs with imag
%     part < 0. It sets anything satisfying this condition to 0. Additionally,
%     anything with imag part between 0 and c has real part set to 0.

global handles;
opt = handles.common.InputSchema({'cut_bias','boundary'}, {false,false}, [], varargin{:});
csqrt = handles.shapelab.common.positive_angle_square_root;

if opt.boundary
  flags = imag(z)<0;
  z(flags) = 0;
  flags = imag(z)<c & abs(real(z))<1e-2;
  %z(flags) = imag(z(flags));
end

tol = 1e-6;
rflags = abs(imag(z))<tol & real(z)<tol;
iflags = (abs(real(z))<tol) & (imag(z)>=0) & (imag(z)<c+tol);

w = csqrt(z.^2 + c.^2,'cut_bias', true);
if not(opt.cut_bias)
  % Preserve x-axis sign + change y-branch
  w(rflags | iflags) = -w(rflags | iflags);
else
  % Just preserve sign on x-axis:
  w(rflags) = -w(rflags);
end
