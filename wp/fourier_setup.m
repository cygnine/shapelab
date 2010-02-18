function[norm_stuff] = fourier_setup(N, varargin)
% fourier_setup -- Precomputations for the WP norm using Fourier Series
%
% norm_stuff = fourier_setup(N, [[filter_coeffs]])
%
%     For some natural number N, this function precomputes some quantities for
%     the `online' computation of the WP norm and computation of Lv using
%     Fourier Series. This assumes that the FFT is to be used, and that the
%     theta interval is [0, 2*pi].
%
%     If filter_coeffs is given, then the filter (which is given as a modal
%     vector of attentuation coefficients) is input.

persistent integer_range ffft_overhead spdiag hilbert_multipliers
if isempty(integer_range)
  from speclab.common import integer_range
  from speclab.fourier.fft import ffft_overhead
  from speclab.fourier import hilbert_multipliers
  from labtools import spdiag
end

if nargin<2
  filter_coeffs = ones([N 1]);
else
  filter_coeffs = varargin{1};
end

norm_stuff.fft_data = ffft_overhead(N, 'shift', pi);
ks = integer_range(N).*filter_coeffs;
norm_stuff.mode_weights = (ks.^3 - ks)/(2*pi);
norm_stuff.mode_weights(ks<2) = 0;
norm_stuff.mode_weights = norm_stuff.mode_weights.';

norm_stuff.lv_weights = spdiag(hilbert_multipliers(N).*((i*ks.^3) - (i*ks))/(2*pi));
