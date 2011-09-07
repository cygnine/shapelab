function[T] = trigonometric_matrix(q, N, d)
% trigonometric_matrix -- The trigonometric interpolation matrix
%
% T = trigonometric_matrix(q, N, [[d=0]])
%
%     Given a vector of Nv positions q, this function returns the
%     Vandermonde-like matrix T that interpolates the basis
%
%     cos(2 x), sin(2 x), ..., cos(N x), sin(N x)
%
%     onto the locations q. The d'th derivatives of the entries are returned. 

if nargin < 3
  d = 0;
end

if N < 2
  error('The input N must be greater than 1');
end

[Ns, q] = meshgrid(2:N, q);

% Doing complex may be faster:
T = [real(exp(i*Ns.*q)) imag(exp(i*Ns.*q))];

if d > 0
  Ns = (2:N).';
  factors = [Ns; Ns].^d;
  T = T*diag(factors);

  if mod(d,2)==0 % no swapping but multiply by sign
    T = T*(-1)^(d/2);
  else  % swap and multiply by sign
    % flip cos and sin
    T = [T(:,N:end) T(:,1:N-1)];

    % multiply cos cols by -1 if d = 1, 5, 9, ...
    if mod((d-1)/2,2)==0
      T(:,1:N-1) = -T(:,1:N-1);
    else % multiply sin cols by -1 if d = 3, 7, 11, ...
      T(:,N:end) = -T(:,N:end);
    end
  end
end
