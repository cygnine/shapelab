function[w] = inverse_base_conformal_map(z,a)
% [w] = inverse_base_conformal_map(z,a)
%
%     Referring to [1], evaluates the inverse of the conformal mapping function
%     f_a, which is a basic building block for the `geodesic algorithm'
%     conformal mapping technique. The complex, scalar parameter is a, and the
%     (array) complex input is z. 
%
%  [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%       mapping", 2006.

global handles;
csqrt = handles.shapelab.common.positive_angle_square_root;
moebius = handles.shapelab.common.moebius;
zipup_at_c = handles.shapelab.common.zipup_at_ic;

assert(length(a)==1, 'Error: not coded for vector-valued parameter a');

% Intermediate points: see [1]
b = abs(a)^2/real(a);
c = abs(a)^2/imag(a);

w = zipup_at_c(z,c,'boundary', true);
w = moebius(w, [1 0; 1/b 1]);
