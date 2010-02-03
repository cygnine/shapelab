function[H] = inverse_map(G)
% inverse_map -- Calculates the inverse of a Moebius map
%
% H = inverse_map(G)
%
%     Given the 2x2 matrix representation of a Moebius map, this function
%     calculates the matrix form of the inverse and returns it.

H = [G(2,2), -G(1,2); -G(2,1), G(1,1)];
