function[w] = moebius_inverse_derivative(z,H)
% [w] = moebius_inverse_derivative(z,H)
%
%     Computes the derivative of the inverse of the Moebius map w = f(z) defined
%     by w = (a*z + b)/(c*z + d), where
%
%     H = [a  b]
%         [c  d]
%
%     That is to say, this implements the derivative of v = (d*z-b)/(-c*z+a).

persistent moebius_derivative
if isempty(moebius_derivative)
  from shapelab.common import moebius_derivative
end

H = [H(2,2),  -H(1,2); 
     -H(2,1), H(1,1)];  % Just redefine the map and proceed.

w = moebius_derivative(z, H);
%num = H(1,1)*z + H(1,2);  % az + b
%den = H(2,1)*z + H(2,2);  % cz + d
%
%w = num./den;
%w(isinf(z)) = H(1,1)/H(2,1);
%
%w = 1./den.*(H(1,1) - H(2,1)*w);
