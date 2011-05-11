function[G] = greens_matrix(q, d)
% greens_matrix -- The Green's function (Gram) matrix for teichons
%
% G = greens_matrix(q, [[d=0]])
%
%     Given a vector of N teichon positions q, this function returns the green's
%     function matrix G with entries:
%
%     G^{(d)}(i,j) = greens_function(q(i) - q(j)),
%
%     where G^{(d)} is the d'th derivative of the WP Dirac Green's function.
%     (see shapelab.wp.greens_function)

persistent greens_function 
if isempty(greens_function)
  from shapelab.wp import greens_function
end

if nargin < 2
  d = 0;
end

[x,y] = meshgrid(q(:), q(:));
G = greens_function((x-y)', d);
