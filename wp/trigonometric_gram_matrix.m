function[TG] = trigonometric_gram_matrix(N)
% trigonometric_gram_matrix -- WP Gram matrix for the trig basis
%
% TG = trigonometric_gram_matrix(N)
%
%     Returns the WP inner product Gram matrix for the basis
%
%     cos(2 x), sin(2 x), ..., cos(N x), sin(N x)
%
%     for N > 1. This is a diagonal matrix with entries
%
%     3, 3, ..., (N^3-N)/2, (N^3-N)/2

if N < 2
  error('The input N must be greater than 1');
end

factors = [repmat((2:N).', [1 2])];
factors = (factors(:).^3 - factors(:))/2;

TG = diag(factors);
