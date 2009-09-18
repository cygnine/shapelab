function[z] = shape_preprocessing(z,varargin)
% shape_preprocessing -- operations to prepare a shape for zipper
%
% [z] = shape_preprocessing(z)
%
%     A shape with vertices z may have the following characteristics, making it
%     unsuitable for zipper:
%       - Nonunique nodes
%       - Intersecting nodes (in this case, nothing is fixed, but a warning is
%         issued)

% Make nodes unique:
[garbage, inds, garbage] = unique(z);
z = z(sort(inds));

ztemp = [z(:); z(1)];
N = length(z);
tol = 1e-8;

% Pairwise search over all line segments for intersection. This is wasteful, but
% *shrug*.
fail = false;
for q = 1:N
  p1 = ztemp(1);
  p2 = ztemp(2);

  for qq = 1:(N-1)
    q1 = ztemp(qq+1);
    q2 = ztemp(qq+2);

    A = [real(q1-q2) real(p2-p1); imag(q1-q2) imag(p2-p1)];
    b = [real(q1-p1); imag(q1-p1)];

    if cond(A)<1e14
      intersect_point = inv(A)*b; s = intersect_point(1); t = intersect_point(2);
      if ((s>tol) & (s<1-tol)) & ((t>tol) & (t<1-tol))
        fail = true;
        break;
      end
    end
  end

  if fail
    warning('The given shape has intersecting segments...the approximation will not work');
    adsf
    break;
  end

  ztemp = circshift(ztemp,-1);
end
