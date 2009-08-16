function[w] = base_conformal_map(z,a,varargin)
% [w] = base_conformal_map(z,a,{cut_bias=true})
%
%     Referring to [1], evaluates the conformal mapping function f_a, which is a
%     basic building block for the `geodesic algorithm' conformal mapping
%     technique. The complex, scalar parameter is a, and the (array) complex
%     input is z. 
%
%     The optional input cut_bias should be toggled depending on how the map is
%     defined. It specifies which side of the csqrt branch the positive x-axis
%     ends up on.
%
%  [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%  mapping", 2006.

global handles;
moebius = handles.shapelab.common.moebius;
unzip_at_c = handles.shapelab.common.z_unzip_at_ic;

assert(length(a)==1, 'Error: not coded for vector-valued parameter a');

% Intermediate points: see [1]
b = abs(a)^2/real(a);
c = abs(a)^2/imag(a);

w = moebius(z, [1 0;-1/b 1]);
w = unzip_at_c(w,c,'boundary',true);
