function[N, N_teeth, z] = assert_enough_points(self, z, varargin)
% assert_enough_points -- Checks for enough points for a shape map
%
% [N, N_teeth, z] = assert_enough_points(z, {other, points, appended, to, z})
%
%     Performs checking to ensure that enough points are given to continue with
%     an unzipping algorithm. For Geodesic-type maps, we need at least 3 points.

z = z(:);
N = length(z);
assert(N>2, 'Error: the geodesic algorithm requires at least three points');
N_teeth = N-2;

% The following does z = [z; varargin{1}; varargin{2}; ...]
z = cat(1, z, varargin{:});
