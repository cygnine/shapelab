function[WP] = fourier_driver(fx, norm_stuff)
% fourier_driver -- Computes the WP norm
%
% WP = fourier_driver(fx, norm_stuff)
%
%     Given equispaced evaluations fx of a function, this computes the WP norm
%     using the fft information in norm_stuff. See fourier_setup.

persistent ffft
if isempty(ffft)
  from speclab.fourier.fft import ffft_online as ffft
end

WP = norm_stuff.mode_weights*abs(ffft(fx, norm_stuff.fft_data)).^2;
