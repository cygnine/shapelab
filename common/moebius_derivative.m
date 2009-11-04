function[w] = moebius_derivative(z,H)
% [w] = moebius_derivative(z,H)
%
%     Computes the derivative of the Moebius map w = f(z) defined by 
%     w = (a*z + b)/(c*z + d), where
%
%     H = [a  b]
%         [c  d]

num = H(1,1)*z + H(1,2);  % az + b
den = H(2,1)*z + H(2,2);  % cz + d

w = num./den;
w(isinf(z)) = H(1,1)/H(2,1);

w = 1./den.*(H(1,1) - w);
