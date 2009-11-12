function[v,w] = oblique_slit_unzip(z,a,varargin)
% [v,w] = oblique_slit_unzip(z,a,{point_id=zeros(size(z))})
%
%     Implements the inverse of oblique_slit_zipup. This function is defined by
%     Marshall in [1]. It is defined as the inverse of the explicit function
%     oblique_slit_zipup.
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

persisten input_schema newton bisection cpow
if isempty(input_schema)
  from labtools import input_schema
  from labtools.rootfind import bisection
  from labtools.rootfind import newton_raphson as newton
  from shapelab.common import positive_angle_exponential as cpow
end

opt = input_schema({'point_id'}, {zeros(size(z))}, [], varargin{:});

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
C = abs(a)/(p^p*q^q);

% The 'interior' case
if any(interior)
  vinfs = not(isfinite(z));
  interior = interior & not(vinfs);
  v(vinfs) = Inf;
  w(vinfs) = Inf;

  % Newton iteration. Initial guess from Marshall:
  v0 = z + 2*p-1 + p*q./(2*z) + ...
       (1-2*p)*p*q./(3*z.^2);

  % Points that are "far": no funny opening-up
  interior_infinity = interior & abs(z)>abs(9/8*a);
  fnorm = z(interior_infinity)/C;
  f = @(x) (x-p./fnorm).^p.*(x+q./fnorm).^q;
  df = @(x) p*((x+q./fnorm)./(x-p./fnorm)).^q + q*((x-p./fnorm)./(x+q./fnorm)).^p;
  %f = @(x) cpow(x-p./fnorm, p).*cpow(x+q./fnorm, q);
  %df = @(x) p*cpow((x+q./fnorm)./(x-p./fnorm), q) + q*cpow((x-p./fnorm)./(x+q./fnorm), p);
  if any(interior_infinity)
    [v(interior_infinity),flag] = newton(v0(interior_infinity)./(fnorm), f, df, ...
      'F', ones(size(z(interior_infinity))),'fx_tol',0,'x_tol',0,'maxiter',30);
    if any(abs(f(v(interior_infinity)) - 1)>1e-8)
      error('Newton''s method didn''t converge');
    end
    v(interior_infinity) = v(interior_infinity).*fnorm;
  end

  % Points that are "near", and "closer to 0"
  interior_0 = interior & (angle(z)/pi < p) & not(interior_infinity);
  % "Open up" this sector
  fnorm = z(interior_0);
  f = @(x) (fnorm/C).^(-1/p).*(x-p).*(x+q).^(q/p);
  df = @(x) (fnorm/C).^(-1/p).*(x+q).^(q/p) + (q/p)*(x-p).*(x+q).^((q-p)/p);
  %f = @(x) (x-p./fnorm).*cpow(x+q./fnorm, q/p);
  %df = @(x) cpow(x+q./fnorm, q/p) + q/p*(x-p./fnorm).*cpow(x+q./fnorm, (q-p)/p);
  if any(interior_0)
    z0 = fnorm + 2*p-1;
    [v(interior_0), flag] = newton(z0, f, df, ...
      'F', ones(size(fnorm)), 'fx_tol', 0, 'x_tol', 0, 'maxiter', 70);
    
    if any(abs(f(v(interior_0)) - 1)>1e-8)
      error('Newton''s method didn''t converge');
    end
    %v(interior_0) = v(interior_0).*fnorm;
  end

  interior_pi = interior & (angle(z)/pi > p) & not(interior_infinity);
  % "Open up" this sector
  fnorm = z(interior_pi);
  %f = @(x) (x-p./fnorm).^(p/q).*(x+q./fnorm);
  %df = @(x) (x-p./fnorm).^(p/q) + p/q*(x-p./fnorm).^((p-q)/q).*(x+q./fnorm);
  f = @(x) (fnorm/C).^(-1/q).*(x-p).^(p/q).*(x+q);
  df = @(x) (fnorm/C).^(-1/q).*((x-p).^(p/q) + (x+q).*(p/q).*(x-p).^(p/q));
  if any(interior_pi)
    %z0 = fnorm + 2*p-1;
    %z0 = v0(interior_pi);
    z0 = fnorm;
    [v(interior_pi), flag] = newton(z0, f, df, ...
      'F', ones(size(fnorm)), 'fx_tol', 0, 'x_tol', 0, 'maxiter', 70);

    if any(abs(f(v(interior_pi)) - 1)>1e-8)
      error('Newton''s method didn''t converge');
    end
  end
end

% If we're on gamma, then the solution on the real line satisfies (x-p)<0 and
% (x+q)>0. We know the phase, so let's just do real-valued iterations. For the
% initial guesses, we'll just use the bounding intervals, [-q,0] and [0,p]
% TODO: can we use bisection?
if any(gamma)
  f = @(x) (p-x).^p.*(x+q).^q;
  %df = @(x) -p*((x+q)./(p-x)).^q + q*((p-x)./(x+q)).^p;

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

  % Newton iteration...fails for a small subset of values
  %v0 = (z(boundary)+ 2*p-1)/C;
  %v0(z(boundary)>0) = max(z(boundary),p+1e-3);
  %v0(z(boundary)<0) = min(z(boundary),-q-1e-3);
  %v(boundary) = newton(v0,f,df,'F',z(boundary)/C);

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

    %[v(boundary),temp] = bisection(v1,v2,f,'F',z(boundary)/C,'x_tol',0,'fx_tol',1e-12);
    % The derivative blows up at z=p...so we can only safely specify x_tol
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
