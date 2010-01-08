function[N, N_teeth, z] = assert_enough_points(z, type, varargin)
% assert_enough_points -- Checks for enough points for a shape map
%
% [N, N_teeth, z] = assert_enough_points(z, type, {other, points, appended, to, z})
%
%     Performs checking to ensure that enough points are given to continue with
%     an unzipping algorithm. The input type is a string that can take the
%     following values:
%
%        'geodesic'
%        'slit'
%        'zipper'
%        'zipper_weld'
%        'loewner'

z = z(:);
N = length(z);

switch lower(type)
case 'geodesic'
  assert(N>2, 'Error: the geodesic algorithm requires at least three points');
  N_teeth = N-2;
case 'slit'
  assert(N>2, 'Error: the slit algorithm requires at least three points');
  N_teeth = N-3;
case 'zipper'
  assert(N>3, 'Error: the zipper algorithm requires at least four points');
  assert(mod(N,2)==0, 'Error: the zipper algorithm requires an even number of points');
  N_teeth = (N-2)/2 - 1;
case 'zipper_weld'
  assert(N>2, 'Error: the zipper-weld algorithm requires at least three points');
  N_teeth = (N-2)/2 - 1;
case 'loewner'
  assert(N>2, 'Error: the loewner algorithm requires at least three points');
  N_teeth = N-3;
otherwise
  error(['Unrecognized algorithm specification: ' type]);
end

% The following does z = [z; varargin{1}; varargin{2}; ...]
z = cat(1, z, varargin{:});
