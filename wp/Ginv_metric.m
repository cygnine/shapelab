function[Gi] = Ginv_metric(theta)
% Ginv_metric -- Computes the inv(G) metric matrix
%
% Ginv = Ginv_metric(theta)
%
%     Returns the inverse of the spd Green's function matrix on the nodes theta.

persistent greens_function
if isempty(greens_function)
  from shapelab.wp import greens_function 
end

N = length(theta);
[thetai, thetaj] = meshgrid(theta, theta);
G = greens_function(thetai-thetaj);

R = chol(G);
Ri = inv(R);
Gi = Ri*Ri';
