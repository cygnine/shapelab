function[k,s] = extend_line_to_circle(k, s, z, newz, varargin)
% extend_line_to_circle -- `Adds' a point to an existing line to make it a circle
%
% [k,s] = extend_line_to_circle(k, s, z, newz, {side=true})
%
%     The scalars k, s, and the complex-valued 2-vector z define a line. (See
%     shapelab.curves.curvature_coefficients, the component z(2) is actually
%     redundant.) This function morphs this representation into that of a circle
%     via specification of a new point newz. The new circle interpolates z and
%     newz, and the outputs are two-vectors k and s specifying the new circle.
%
%     The optional boolean input side determines is the order of interpolation
%     is [z newz] or [newz z]. The default (side=true) is the former.

persistent circle_coefficients strict_inputs
if isempty(circle_coefficients)
  from shapelab.curves import circle_coefficients
  from labtools import strict_inputs
end

% Let's try to get around using strict_inputs for speed
if nargin==4
  side = true;
else 
  side = varargin{2};
end

assert(length(s)==1, 'Error: s must be a scalar');
assert(length(k)==1, 'Error: k must be a scalar');

% Don't actually use z(2)
z2 = z(1) + s*exp(i*k);

% Now z(1), z2 and newz are three points. Find the circle:
if side
  z = [z(1); z2; newz];
else
  z = [newz; z(1); z2];
end

[k,s] = circle_coefficients(z);
