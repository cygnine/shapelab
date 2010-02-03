function[norm_stuff] = fourier_setup(N)
% fourier_setup -- Precomputations for the WP norm using Fourier Series
%
% norm_stuff = fourier_setup(N)
%
%     For some natural number N, this function precomputes some quantities for
%     the `online' computation of the WP norm using Fourier Series. This assumes
%     that the FFT is to be used, and that the theta interval is [0, 2*pi].

persistent integer_range ffft_overhead spdiag
if isempty(integer_range)
  from speclab.common import integer_range
  from speclab.fourier.fft import ffft_overhead
  from labtools import spdiag
end

norm_stuff.fft_data = ffft_overhead(N, 'shift', pi);
ks = integer_range(N);
norm_stuff.mode_weights = (ks.^3 - ks)/(2*pi);
norm_stuff.mode_weights(ks<2) = 0;
norm_stuff.mode_weights = norm_stuff.mode_weights.';
