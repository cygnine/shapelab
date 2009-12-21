function[ind] = indexing(scale)
% indexing -- Anonymous function handles for scaling
%
% ind = indexing(scale)
%
%     Helper function for point-value differences. Based on input `scale',
%     returns an anonymous function `ind', a function of one positive integer k
%     that determines how to index `scale k'. 

switch scale
case 'dyadic'
  ind = @(k) 2^k;
case 'linear'
  ind = @(k) k+1;
end
