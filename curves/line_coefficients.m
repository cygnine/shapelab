function[k,s] = line_coefficients(z)
% line_coefficients -- Determines curvature coefficients for a line
%
% [k,s] = line_coefficients(z)
%
%     Given a 2-vector with complex entries z, this function determines the
%     scalars k and s such that p(t) = k is a primitive of the 0 curvature and
%     satisfies
%
%     z(2) = z(1) + \int_0^s(1) exp(i*p(t)) dt

s = abs(z(2)-z(1));
k = angle(z(2)-z(1));
