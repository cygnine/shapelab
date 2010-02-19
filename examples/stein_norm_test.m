% Script to test approximation via the Stein formulas of the WP metric

clear 
close all

from shapelab.wp import stein_norm

% The first run of stein_norm(.) always takes more time than successive runs.
% Some things must be precomputed and stored. If you want to increase the
% accuracy of the method, you can mess around with the parameter N given in
% stein_norm

qs = 1:20;
for q = qs
  v = @(x) exp(i*q*x);
  E(q) = stein_norm(v);
end

% exact is n^3-n
fprintf('The computed WP norms of complex exponentials (left). Exact values are on the right.\n');
[E.' (qs.^3 - qs).']
