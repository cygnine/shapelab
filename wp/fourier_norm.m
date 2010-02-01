function[WP] = fourier_norm(fx)
% fourier_norm -- Computes the WP norm using equispaced samples of a function
%
% WP = fourier_norm(fx)
%
%     Uses an equispaced Fourier quadrature rule to compute the WP norm of a
%     function. The input fx is a length-N vector containing the function
%     evaluations -- it is assumed that the nodal locations are the equispaced
%     Fourier quadrature rule output from speclab.fourier.quad.gauss_quadrature.

persistent fourier_setup fourier_driver
if isempty(fourier_setup)
  from shapelab.wp import fourier_setup fourier_driver
end

norm_stuff = fourier_setup(length(fx));
WP = fourier_driver(fx, norm_stuff);
