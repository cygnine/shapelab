function[norms] = fourier_coeffs_norm(coeffs)
% fourier_coeffs_norm - Calculates the WP norm given Fourier coefficients
%
% wp_norm = fourier_coeffs_norm(coeffs)
%
%     Given Fourier coefficients of a function expansion, this computes the
%     Weil-Petersson norm of that expansion.  When coeffs is a vector:
%     L_coeffs_k <---- |k^3 - k| coeffs_k, and then the norm is computed as
%     sum(L_coeffs_k.*coeffs_k). The indexing is assumed to be negatively-biased
%     integer indexing (as in speclab.common.integer_range).
%
%     When coefffs is a matrix, the vector operation is performed on each
%     column. 

persistent integer_range 
if isempty(integer_range)
  from speclab.common import integer_range
end

[m,n] = size(coeffs);
if (m==1) | (n==1)
  coeffs = coeffs(:);
  N = max([m,n]);
else
  N = m;
end

ks = integer_range(N);

factors = abs(ks.^3 - ks);
norms = 1/(2*pi)*factors.'*abs(coeffs).^2;

% The 1/(2*pi) comes in because I assume the coefficients are for basis
% functions 1/sqrt(2*pi)*exp(i*k*theta), whereas the k^3-k coefficient
% dependence is for basis functions exp(i*k*theta).
