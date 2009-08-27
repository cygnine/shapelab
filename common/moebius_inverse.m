function[v] = moebius_inverse(z,H)
% [v] = moebius_inverse(z,H)
%
%     Implements the *inverse* of the Moebius map w = f(z) defined by
%     w = (a*z + b)/(c*z+d), where
%
%     H = [a  b]
%         [c  d]
%
%     That is to say, this implements v = (d*z-b)/(-c*z+a).

v = (H(2,2)*z - H(1,2))./(-H(2,1)*z + H(1,1));
v(isinf(z)) = -H(2,2)/H(2,1);
