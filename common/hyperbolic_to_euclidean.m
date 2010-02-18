function[w] = hyperbolic_to_euclidean(z)
% hyperbolic_to_euclidean -- Changes to euclidean distance from 0
%
% w = euclidean_to_hyperbolic(z)
%
%     If the input z has the form rho*exp(i*theta), where rho is the hyperbolic
%     distance from 0, the output w has the same angle, but the magnitude is now
%     the euclidean distance from 0.

rho = abs(z);
w = tanh(rho).*exp(i*angle(z));
