function[w] = normal_linear_slit_zip(a, z)
% normal_linear_slit_zip -- Zips up a linear slit normal to the real axis
%
% [w, w_slit_left, w_slit_right] = normal_linear_slit_zip(a, z)
%
%     The exact (inverse) solution to the Loewner equation for a linear slit
%     that is perpendicular to the axis at the location x = real(a).
%
%                                                                a
%                                                               x
%                                                               |                    
%          z-plane                                     w-plane  |                
%                                         ---->                 |                
%                                                               |
%   real(a)-imag(a)    real(a)                                  |                
%         _________o-----x-----o________        ________________o________________
%
%
%     Unlike this function's inverse (normal_linear_slit_unzip), there is no
%     ambiguity regarding which branch points are mapped to.

persistent cexp
if isempty(cexp)
  from shapelab.common import positive_angle_exponential as cexp
end

w = zeros(size(z));
ra = real(a);
ia = imag(a);

% There is one exception that must be made: any values that are real and less
% than real(a) - imag(a)
reals = (imag(z)==0);
critical_negative = (real(z) < (ra-ia));
flags = reals & critical_negative;

w = cexp(((z - ra).^2 - ia^2), 1/2) + ra;

% Here's the exception: it comes about because cexp can't differentiate between
% polar theta = 0 and theta = 2*pi.
w(flags) = ra - (w(flags) - ra);
