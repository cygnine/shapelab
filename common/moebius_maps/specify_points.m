function[H] = specify_points(points, images)
% specify_points -- Determines Moebius maps from point-image pairs
%
% H = specify_points(points, images)
%
%     Both inputs `points' and `images' are 3-vectors. The output Moebius map H
%     takes points(i) to images(i) for i = 1, 2, 3.

persistent map_inf_to preprocess inverse_map
if isempty(map_inf_to)
  from shapelab.common.moebius_maps import map_infinity_to as map_inf_to
  from shapelab.common.moebius_maps import finitize_points as preprocess
  from shapelab.common.moebius_maps import inverse_map
end

points = points(:);
images = images(:);

% Get rid of Infs:
[points, pre_map] = preprocess(points);
[images, post_map] = preprocess(images);
post_map = inverse_map(post_map);

% First try to solve for:
% w = (a*z + b)/(z + d)      (1)
A = [points(:), ones([3 1]), -images(:)];

tol = 1e-20;
if abs(det(A))<tol % Then guess (1) was wrong.
  % The guess a new form:
  % w = (a*z+b)

  % If all three `points' are the same:
  if all(abs(diff([points(:); points(1)]))<tol)
    if not(all(abs(diff([images(:); images(1)])))<tol)
      fprintf('The given points/images don''t have a well-posed Moebius map\n');
    else
      H = [0, images(1); 0, 1];  % The trivial (constant) map
    end
  else % Proceed
    % Just find two of the points that aren't coincident
    temp = diff(points);
    if abs(temp(1))>tol
      inds = [1; 2];
    elseif abs(temp(2))>tol
      inds = [2; 3];
    else
      inds = [1; 3];
    end

    A = [points(inds), ones([2 1])];
    coeffs = inv(A)*images(inds);

    H = [coeffs(1), coeffs(2); 0, 1];
  end

else  % Then guess (1) is right -- proceed to invert system
  coeffs = inv(A)*(points(:).*images(:));
  H = [coeffs(1) coeffs(2); 1 coeffs(3)];
end

H = post_map*H*pre_map;
