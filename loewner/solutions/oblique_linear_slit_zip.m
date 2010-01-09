function[w] = oblique_linear_slit_zip(a, alpha, z, varargin)
% oblique_linear_slit_zip -- Zips up a linear slit that is oblique to the real axis
%
% [w] = oblique_linear_slit_unzip(a, alpha, z, {{ Cn=1 }});
%
%     The exact solution to the Loewner equation for a linear slit that is
%     oblique to the axis at the location x = c (derived from alpha and a),
%     where alpha is the angle (radians) made with the axis.
%
%                                                                   a
%                                                                    x
%                                                                   /
%         z-plane                               w-plane            /            
%                                     ----->                      /             
%                                                                /  
%                        c                                   c  / ) alpha       
%        _________o-----x---o________          ________________o__)_____________
%
%
%     This map is the inverse of oblique_linear_slit_unzip.
%
%     Cn is a normalization constant -- it should be consistent when called from
%     oblique_linear_slit_unzip

if nargin>3
  Cn = varargin{1};
else
  Cn = 1;
end

w = zeros(size(z));
ra = real(a);  % The point on the real axis to unzip to 
ia = imag(a);
c = a - ia/sin(alpha)*exp(i*alpha);

% Marshall's terminology:
p = alpha/pi;
q = 1-p;
C = abs(a-c)/(p^p*q^q);

f = @(x) C*((x-c)/Cn/C-p).^p.*((x-c)/Cn/C+q).^q;

w = f(z) + c;

% Once again, sometimes we get machine eps crap:
flags = imag(w)<0;
w(flags) = real(w(flags));

w(isinf(z)) = Inf;
