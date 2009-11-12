function[w] = inverse_base_conformal_map(z,c,a,varargin)
% [w] = inverse_base_conformal_map(z,c,a,{point_id=zeros(size(z)),cut_magnitude=abs(a))
%
%     Referring to [1], evaluates the inverse of the conformal mapping function
%     f_a, which is a basic building block for the `zipper algorithm'
%     conformal mapping technique. The complex, scalar parameters are a and c, and the
%     (array) complex input is z. The parameter a is the 'tip' of the circular
%     arc, whereas c is between 0 and a on the arc.
%
%     The optional input point_id has the same size as z and each entry takes
%     three possible values:
%     0: The point is some point in \mathbb{H}\backslash\mathbb{R}.
%     1: The point is located on \mathbb{R}
%     2: The point is located on \mathbb{R}, inside [p-1,p]. (See code for
%        definition of p.)
%
%  [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%       mapping", 2006.

persistent zipup moebius moebius_inv input_schema
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import moebius
  from shapelab.common import moebius_inverse as moebius_inv
  from shapelab.common import slit_zipup_to_a as zipup
end
opt = input_schema({'point_id'}, ...
      {zeros(size(z),'int8')}, [],varargin{:});

assert(length(a)==1, 'Error: not coded for vector-valued parameter a');
assert(length(c)==1, 'Error: not coded for vector-valued parameter c');

delta_real = real(a) - real(c);
if abs(delta_real)<1e-14
  d = a;  % Straight line segment...we *could* actually use geodesic here since
          % it's the same conformal map
  m = eye(2);
else
  b = (imag(c)*abs(a)^2 - imag(a)*abs(c)^2)/(real(a)*imag(c) - real(c)*imag(a));
  m = [1 0; -1/b 1];
  d = moebius(a,m);
end

w = zipup(z,d,'point_id',opt.point_id);
w = moebius_inv(w,m);
