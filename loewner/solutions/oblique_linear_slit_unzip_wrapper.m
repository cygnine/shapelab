function[z, z_interior, z_exterior] = oblique_linear_slit_unzip_wrapper(a, alpha, z, z_interior, z_exterior)
% oblique_linear_slit_unzip_wrapper -- Syntactic sugar for oblique_linear_slit_unzip
% 
% [z, z_interior, z_exterior] = ...
%   oblique_linear_slit_unzip_wrapper(a, z, z_interior, z_exterior)
%
%     A wrapper for zipper-like functions.

persistent unzip
if isempty(unzip)
  from shapelab.loewner.solutions import oblique_linear_slit_unzip as unzip
end

Next = length(z_exterior(:));
Nint = length(z_interior(:));
z_slit = [z_interior(:); z_exterior(:)];

[z, temp2, temp3] = unzip(a, alpha, z, z_slit);

if Nint>0
  z_interior = temp2(1:Nint);
end
if Next>0
  z_exterior = temp3((Nint+1):end);
end
