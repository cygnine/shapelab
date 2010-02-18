function[z] = D_to_H_cover(z)
% D_to_H_cover -- The inverse of the covering exp(i*z) from H to D
%
% z = D_to_H_cover(z)
%
%     Implements the inverse of the covering z -----> exp(i*z)

z = -i*log(z);
