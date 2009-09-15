function[v,w] = slit_unzip_from_a(z,a,varargin)
% [v,w] = slit_unzip_from_a(z,a,{point_id=zeros(size(z))})
%
%     Implements the inverse of slit_zipup_to_a. This function is defined by
%     Marshall in [1]. It is defined as the inverse of the explicit function
%     slit_zipup_to_a.
%
%     This function computes the inverse of the map
%          z(w) = C*(w-p)^p*(w+q)^q,
%     where p = angle(a)/pi, q = 1-p, and C = |a|/(p^p*q^q). Without loss we
%     normalize a to have norm 1.
%
%     The outputs v and w are identical except for any z on gamma (defined
%     below), when they represent the negative and positive solutions,
%     respectively. point_id takes the following values:
%
%     0: The point is somewhere in \mathbb{H}\backslash{\gamma}, where \gamma is
%        the line segment (0,0) --a
%     1: The point is located on gamma
%     2: The point is located on \mathbb{R}
%
% [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%      mapping", 2006.
%
%     TODO: currently using bisection on 2...get Newton's method to converge

global handles;
opt = handles.common.input_schema({'point_id'}, {zeros(size(z))}, [], varargin{:});
newton = handles.common.rootfind.newton_raphson;
bisection = handles.common.rootfind.bisection;
cpow = handles.shapelab.common.positive_angle_exponential;

interior = opt.point_id==0;
gamma = opt.point_id==1;
boundary = opt.point_id==2;

v = zeros(size(z));
w = zeros(size(z));

% pre-processing
p = angle(a)/pi;
q = 1-p;
if p<=0 | p>=1
  error('The input a must be in the upper half plane');
end

% Normalize a (and the inputs):
z = z/abs(a);
% Now abs(a)=1:
C = 1/(p^p*q^q);

% The 'interior' case
if any(interior)
  % First: exempt infinity from any manipulations
  vinfs = not(isfinite(z));
  interior = interior & not(vinfs);
  v(vinfs) = Inf;
  w(vinfs) = Inf;

  % On the interior, do Newton's method, but take initial guess as following:
  % for angle(z)<pi*p, use z^(0.5/p) to map z to upper-right quadrant. Then use
  % the solution for symmetric_unzip_from_ic, and then raise to power z^(p/0.5).
  % Do a similar thing for angle(z)>pi*p.
  small_angles = interior & angle(z)/pi<p;
  large_angles = interior & not(small_angles);

  % Initial guess for small angles:
  v0 = zeros(size(z));
  z(small_angles) = z(small_angles).^(0.5/p);
  z(large_angles) = ((z(large_angles)*exp(-i*pi*p)).^(0.5/q)).*exp(i*pi/2);
  v0(interior) = slit.symmetric_unzip_from_ic(z,1);
  v0(small_angles) = v0(small_angles).^(p/0.5);
  v0(large_angle) = ((v0(large_angles).*exp(-i*pi/2)).^(q/0.5)).*exp(i*pi*p);

  % Newton iteration. Initial guess from Marshall:
  %v0 = z + 2*p-1 + p*q./(2*z) + ...
  %     (1-2*p)*p*q./(3*z.^2);

  % Points that are "far": no funny opening-up
  % If the default doesn't work, maybe go the Marshall way and divide by fnorm
  % so that f(w) = 1?
  %fnorm = z(interior)/C;
  fnorm = ones(size(z(interior)));
  f = @(x) C*(x-p).^p.*(x+q).^q;
  df = @(x) C*p*(x-p).^(p-1).*(x+q).^q + C*q*(x-p).^p.*(x+q).^(q-1);
  %f = @(x) (x-p./fnorm).^p.*(x+q./fnorm).^q;
  %df = @(x) p*((x+q./fnorm)./(x-p./fnorm)).^q + q*((x-p./fnorm)./(x+q./fnorm)).^p;
  if any(interior)
    [v(interior),flag] = newton(v0(interior), f, df, ...
      'F', z(interior),'fx_tol',0,'x_tol',0,'maxiter',30);
    if any(abs(f(v(interior)) - z(interior))>1e-8)
      error('Newton''s method didn''t converge');
    end
    %v(interior) = v(interior).*fnorm;
  end
end

% If we're on gamma, then the solution on the real line satisfies (x-p)<0 and
% (x+q)>0. We know the phase, so let's just do real-valued iterations. For the
% initial guesses, we'll just use the bounding intervals, [-q,0] and [0,p]
% TODO: can we use bisection?
if any(gamma)
  f = @(x) (p-x).^p.*(x+q).^q;

  % Negative side
  v0 = -q*ones(size(z(gamma)));
  v1 = 0*v0;
  v(gamma) = bisection(v0,v1,f,'F',abs(z(gamma))/C);

  % positive side
  w0 = 0*v0;
  w1 = p*ones(size(w0));
  w(gamma) = bisection(w0,w1,f,'F',abs(z(gamma))/C);
end

% On the real line -- god this is messy
if any(boundary)

  vinfs = not(isfinite(z));
  boundary = boundary & not(vinfs);
  v(vinfs) = Inf;
  w(vinfs) = Inf;

  v0s = (abs(z)<1e-13);
  boundary = boundary & not(v0s);
  v(v0s) = -q;
  w(v0s) = p;

  % Trust me on this
  f = @(x) sign(x).*(abs(x-p).^p.*abs(x+q).^q);
  %df = @(x) p*abs((x+q)./(x-p)).^q + q*abs((x-p)./(x+q)).^p;

  % Bisection is nearly foolproof
  if any(boundary)
    v1 = zeros(size(z(boundary)));
    v2 = v1;
    flags = z(boundary)>0;
    temp = z(boundary);
    if any(flags);
      v1(flags) = p*ones(size(temp(flags)));
      v2(flags) = (temp(flags)+2*p-1)/C + 1e2;
    end

    if any(~flags)
      v1(~flags) = (temp(~flags)+2*p-1)/C - 1e2;
      v2(~flags) = -q*ones(size(temp(~flags)));
    end

    [v(boundary),temp] = bisection(v1,v2,f,'F',z(boundary)/C,'x_tol',1e-14);
    if temp==1
      error('Bisection didn''t converge');
    end
    w(boundary) = v(boundary);
  end

end

if any(imag(v)<0)
  error('Stuff''s in the lower half-plane');
end
