function[a] = compute_a(g,lambda)
% compute_a -- computes a(g,lambda)
%
% a = compute_a(g)
%
%     Computes the function a(g) defined as 
%
%       a = 1/2*(g.^2) - lambda*g
%
%     This map preserves infinity.

flags = isinf(g);

a = 1/2*g.^2 - lambda*g;
a(flags) = Inf;
