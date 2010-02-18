function[w] = euclidean_to_hyperbolic(z)
% euclidean_to_hyperbolic -- Changes to hyperbolic distance from 0
%
% w = euclidean_to_hyperbolic(z)
%
%     If the input z has the Euclidean form r*exp(i*theta), the output w has the
%     same angle, but the magnitude is now the hyperbolic distance from 0.

r = abs(z);
w = exp(i*angle(z))*1/2.*log((1+r)./(1-r));
