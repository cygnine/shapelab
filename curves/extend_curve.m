function[k,s] = extend_curve(k, s, z, newz, varargin)
% extend_curve -- Extends a curvature interpolation by one point
%
% [k,s] = extend_curve(k, s, z, newz, {side=true})
%
%     Extends a curve defined by the monomial curvature coefficients k at
%     arclengths s by the point newz. The current z values are z, and if 'side' is
%     true the new data is [z; newz], otherwise it is [newz; z]. This ordering
%     is important. 

persistent strict_inputs eval_curvature_samples monomial_shift ...
           extend_line_to_circle newton_solve polar_interpolate
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.curves import eval_curvature_samples extend_line_to_circle ...
                              newton_solve polar_interpolate
  from speclab.monomials import monomial_shift
end

opt = strict_inputs({'side'}, {true}, [], varargin{:});

z = z(:);
N_small_steps = 20;
N = length(k);

if N==1  % Extending a line to a circle is easy:
  [k,s] = extend_line_to_circle(k, s, z, newz, 'side', opt.side);
  return
end

if opt.side  % append
  x = [k(:); 0; s(:); s(end)+abs(newz-z(end))];
  zend = z(1) + eval_curvature_samples(k(:), x(end));
  %z_ends = linspace(zend, newz, N_small_steps+1);
  z_ends = z(end) + polar_interpolate(zend-z(end), newz-z(end), N_small_steps+1);

  for q = 1:(N_small_steps+1)
    ztemp = [z; z_ends(q)];

    [x, converged] = newton_solve(ztemp, x);
    if not(converged)
      k = NaN;
      s = NaN;
      return
    end
  end

  k = x(1:(N+1));
  s = x((N+2):end);

else % prepend
  shift = abs(newz-z(1));
  z1 = z(1) + eval_curvature_samples(k(:),-shift);

  k = monomial_shift(k, shift);

  x = [k(:); 0; shift; s(:)+shift];
  z_ends = linspace(z1, newz, N_small_steps+1);
  z_ends = z(1) + polar_interpolate(z1 - z(1), newz - z(1), N_small_steps+1);

  for q = 1:(N_small_steps+1)
    ztemp = [z_ends(q); z];

    [x, converged] = newton_solve(ztemp, x);
    if not(converged)
      k = NaN;
      s = NaN;
      return
    end
  end

  k = x(1:(N+1));
  s = x((N+2):end);
end
