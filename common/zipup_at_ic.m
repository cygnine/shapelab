function[w] = zipup_at_ic(z,c,varargin)
% [w] = z_unzip_at_ic(z,c,{cut_bias=false,boundary=false})
%
%     Implements the inverse of z_unzip_at_c. I.e., the function sqrt(z^2 -
%     c^2). 

global handles;
tol = 1e-6;
csqrt = handles.shapelab.common.positive_angle_square_root;
opt = handles.common.InputSchema({'boundary'}, {false}, [], varargin{:});

if opt.boundary
  flags = abs(imag(z))<=tol;
  z(flags) = real(z(flags));
end

w = csqrt(z.^2 - c^2,'cut_bias', true);

rflags = abs(imag(z))<tol & real(z)<(-c+tol);
w(rflags) = -w(rflags);
