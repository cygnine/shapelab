function[r,x0] = interpolate_circle(a,b,c)
% interpolate_circle -- interpolates points to form a circle
%
% [r,x0] = interpolate_circle(a,b,c)
%
%     Given three 2-tuples a, b, and c, this function returns Euclidean
%     parameters r and x0 defining any point x on the circle:
%
%       (x(1) - x0(1))^2 + (y(1) - y0(1))^2 = r^2

denom = (a(1)-b(1))*(b(2)-c(2)) - (a(2)-b(2))*(b(1)-c(1));

if abs(denom)<1e-14
  r = Inf;
  x0 = a;
  return
end

am = sqrt(sum(abs(a).^2));
bm = sqrt(sum(abs(b).^2));
cm = sqrt(sum(abs(c).^2));

A = 2*[b(1)-a(1) b(2)-a(2);
       c(1)-b(1) c(2)-b(2)];

rhs = [bm-am; cm-bm];

x0 = inv(A)*rhs;

%x0(1) = bm*(c(1)-a(1)) + am*(b(1)-c(1)) + cm*(a(1)-b(1));
%x0(2) = bm*(c(2)-a(2)) + am*(b(2)-c(2)) + cm*(a(2)-b(2));

%x0 = 1/4*x0.*([1./denom, -1./denom]);

r = sqrt(sum(abs(a(:)-x0(:)).^2));
