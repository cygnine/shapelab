function[m,x0] = interpolation_line(a,b)
% interpolate_line -- interpolates points to form line
%
% [m,x0] = interpolate_line(a,b)
%
%     Given two 2-tuples a and b, the function returns the Euclidean parameters
%     m and x defining any point x on the line via point-slope form:
%
%        x(2) - x0(2) = m*(x(1) - x0(1))

if a(1)==b(1)
  m = Inf;
else
  m = (b(2)-a(2))/(b(1)-a(1));
end

x0 = a;
