function[w,w_slit_left,w_slit_right] = oblique_linear_slit_unzip(a, alpha, z, varargin)
% oblique_linear_slit_unzip -- Unzips a linear slit that is oblique to the real axis
%
% [w, w_slit_left, w_slit_right] = oblique_linear_slit_unzip(a, alpha, z, ...
%                                   {{z_slit = [], Cn=1}})
%
%     The exact solution to the Loewner equation for a linear slit that is
%     oblique to the axis at the location x = c (derived from alpha and a),
%     where alpha is the angle (radians) made with the axis.
%
%                          a
%                           x
%                          /
%      z-plane            /                        w-plane
%                        /                ----->     
%                       /  
%                   c  / ) alpha                                  c            
%     ________________o__)_____________           _________o-----x---o________
%
%
%     The input vector z contains all points that are not on the slit; the map
%     in this case is simply implemented as there is no ambiguity. The optional
%     input z_slit contains points on the slit that are to be unzipped. No
%     testing is performed on z_slit to make sure it actually does lie on the
%     slit. If z_slit is given, then both the `left' and `right' images are
%     given in output.
%
%     If you just want to evaluate the map for points on the slit, set the
%     mandatory input z to the empty array ([]).
%
%     Cn is a normalization constant -- it multiplies the output w by the
%     constant factor Cn.

persistent compute_driving drive newton bisection
if isempty(drive)
  from shapelab.loewner import compute_driving_function as compute_driving
  from shapelab.loewner import pointwise_driving as drive
  from labtools.rootfind import newton_raphson as newton
  from labtools.rootfind import bisection
end

if nargin>3
  z_slit = varargin{1};
  slit = true;
  if nargin>4
    Cn = varargin{2};
  else
    Cn = 1;
  end
else
  slit = false;
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

% First compute the driving function
M = 150; % Number of Forward-Euler steps to take
[s, lambda, dlambda] = compute_driving(linspace(0, a-c, M+1));

% All the regular points:
if not(isempty(z))
  w_guess = drive(s, lambda, dlambda, z-c) - lambda(end) + c;
  % Use the above as an initial guess

  f = @(x) C*(x-p).^p.*(x+q).^q;
  df = @(x) C*p*(x-p).^(p-1).*(x+q).^q + C*q*(x-p).^p.*(x+q).^(q-1);
  w = newton((w_guess-c)/C, f, df, 'F', z-c, 'fx_tol', 1e-12, 'x_tol', 0, 'maxiter', 1e2);

  troubles = find(abs(f(w)-(z-c))>1e-8);
  w(troubles) = newton((w_guess(troubles)-c)/C, f, df, 'F', z(troubles)-c, 'fx_tol', 1e-12, 'x_tol', 0, 'maxiter', 1e2, 'tiptoe', 0.1);
end

% Sometimes imaginary parts get set to -eps
flags = imag(w)<0;
w(flags) = real(w(flags));

% The mapping on the slit can be done via bisection
f = @(x) (p-x).^p.*(x+q).^q;

% In particular, we have real-valued iterations
if slit
  % Negative side
  v0 = -q*ones(size(z_slit));
  v1 = 0*v0;
  w_slit_left = bisection(v0,v1,f,'F',abs(z_slit)/C);

  % positive side
  w0 = 0*v0;
  w1 = p*ones(size(w0));
  w_slit_right = bisection(w0,w1,f,'F',abs(z_slit)/C);
else
  w_slit_left = [];
  w_slit_right = [];
end

% Normalization
w = w*(C*Cn);
w_slit_left = w_slit_left*(C*Cn);
w_slit_right = w_slit_right*(C*Cn);
