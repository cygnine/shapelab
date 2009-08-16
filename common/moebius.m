function[w] = moebius(z,H)
% [w] = moebius(z,H)
%
%     Implements the Moebius map w = f(z) defined by
%     w = (a*z + b)/(c*z+d), where
%
%     H = [a  b]
%         [c  d]

w = (H(1,1)*z + H(1,2))./(H(2,1)*z + H(2,2));
w(isinf(z)) = H(1,1)/H(2,1);
