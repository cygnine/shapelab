function[wp] = teichon_wp(a, b)
% teichon_wp -- Computes the WP norm of a teichon velocity field
%
% wp = teichon_wp(a, b)
%
%     Given a velocity field of the form
%
%       v(theta) = sum_j a_j G(theta-b_j),
%
%     where G is the WP-Green's function, this function computes and returns the
%     WP-norm.

persistent greens_function
if isempty(greens_function)
  from shapelab.wp import greens_function
end

[temp1, temp2] = meshgrid(b,b);
wp = sqrt(a'*greens_function(temp1-temp2)*a);
