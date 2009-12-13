function[ks, ss, z0] = geno_interpolant(z,varargin)
% geno_interpolant -- returns a curve interpolant constructed using GENO
%
% [ks,ss, z0] = geno_interpolant(z,{closed=false, k=2})
%
%     Interpolates the ordered set of the complex-valued points z with a
%     curve that has piecewise polynomial curvature of order (k-2). I.e. it
%     interpolates with circles when k=2, and with straight lines when k=1. If
%     the input closed is set to true then the interpolation is 'periodic' in the
%     sense that z(1) and z(end) are connected via interpolation as well.
%
%     When closed is false: If length(z) is N, the output ks is a k x (N-1)
%     matrix, where column q represents a polynomial of the primitive curvature
%     for segment z(q) -- z(q+1). ss is a 2 x (N-1) matrix where the segment
%     lies between ss(q,1) and ss(q,2). The offset is z0(q).
%
%     When closed is true: If length(z) is N, the output ks is a k x N matrix,
%     where column q=1, 2,...,N-1 represents a polynomial of the primitive
%     curvature for segment z(q) -- z(q+1). Column q=N is for the segment
%     z(N)--z(1). ss is a 2 x (N-1) matrix where the segment lies between
%     ss(q,1) and ss(q,2). The offset is z0(q).

persistent strict_inputs curv_coeffs extend_curve index_mod
if isempty(strict_inputs)
  from labtools import strict_inputs index_mod
  from shapelab.curves import curvature_coefficients as curv_coeffs
  from shapelab.curves import extend_curve
end

inputs = {'closed', 'k'};
defaults = {false, 2};
opt = strict_inputs(inputs, defaults, [], varargin{:});

N = length(z);
choose_left = false;

if not(opt.closed)  % No periodicity
  ks = zeros([opt.k N-1]);
  ss = zeros([2 N-1]);
  z0 = zeros([1 N-1]);

  for q = 1:(N-1);
    % Sucessively build up each connecting segment
    current_stencil = [q; q+1];

    [k,s] = curv_coeffs(z(current_stencil));

    for qq = 2:opt.k;
      left_stencil = [current_stencil(1)-1; current_stencil];
      right_stencil = [current_stencil; current_stencil(end)+1];

      if left_stencil(1)>0
        [k_l,s_l] = extend_curve(k, s, z(current_stencil), z(left_stencil(1)), 'side', false);
      else
        k_l = [];
      end

      if right_stencil(end)<N
        [k_r,s_r] = extend_curve(k, s, z(current_stencil), z(right_stencil(end)), 'side', true);
      else
        k_r = [];
      end

      choose_side();
      assign_ks();

    end

    s = [0; s(:)];
    ks(:,q) = k;
    ss(1,q) = s(current_stencil==q);
    ss(2,q) = s(current_stencil==(q+1));
    z0(q) = z(current_stencil(1));

  end

else  % closed curve

  ks = zeros([opt.k N]);
  ss = zeros([2 N]);
  z0 = zeros([1 N]);

  for q = 1:N; % for each segment
    % Sucessively build up each connecting segment
    current_stencil = index_mod([q; q+1],N);

    [k,s] = curv_coeffs(z(current_stencil));

    break_interpolation = false;

    for qq = 2:opt.k;
      left_stencil = index_mod([current_stencil(1)-1; current_stencil],N);
      right_stencil = index_mod([current_stencil; current_stencil(end)+1],N);

      [k_l,s_l] = extend_curve(k, s, z(current_stencil), z(left_stencil(1)), 'side', false);
      [k_r,s_r] = extend_curve(k, s, z(current_stencil), z(right_stencil(end)), 'side', true);

      if abs(k_l(end))<abs(k_r(end))   % Choose k_l
        choose_left = true;
      else
        choose_left = false;
      end

      choose_side();
      if break_interpolation
        % Pad vectors with zeros and quit
        k = [k; zeros([opt.k-qq+1 1])];
        s = [s; zeros([opt.k-qq+1 1])];
        break
      end
      assign_ks();

    end

    s = [0; s(:)];
    ks(:,q) = k;
    ss(1,q) = s(current_stencil==q);
    ss(2,q) = s(current_stencil==(index_mod(q+1,N)));
    z0(q) = z(current_stencil(1));

  end

end

%%%%%%%%%%%%%%% Nested functions %%%%%%%%%%%%%%% 

function[] = choose_side()
% choose_side -- nested function
%
% choose_side
%
%     Based on values of k_l and k_r, determines which boolean value to set
%     choose_left to.

  if isempty(k_l)
    choose_left = false;
  elseif isempty(k_r)
    choose_left = true;
  elseif isnan(k_l(1))
    choose_left = false;
    if isnan(k_r(1))
      disp('Couldn''t interpolate to either side');
      break_interpolation = true;
    end
  elseif isnan(k_r(1))
    choose_left = true;
  elseif abs(k_l(end))<abs(k_r(end)) % Choose k_l
    choose_left = true;
  else
    choose_left = false;
  end

end

function[] = assign_ks()
% assign_ks -- nested function
%
% assign_ks()
%
%     Just assigns variables k, s, and current_stencil based on the value of
%     choose_left.

  if choose_left
    current_stencil = left_stencil;
    k = k_l;
    s = s_l;
  else
    current_stencil = right_stencil;
    k = k_r;
    s = s_r;
  end

end

end
