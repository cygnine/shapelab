function[points, map] = finitize_points(points)
% finitize_points -- Makes points finite via a Moebius map
%
% [points, map] = finitize_points(points)
%
%     This function makes a collection of points finite via a Moebius map that
%     is a rotation between infinity and some selected point. The output
%     points is the collection of points that have been rotated to be finite,
%     and the matrix Moebius map is the second output.

persistent moebius
if isempty(moebius)
  from shapelab.common import moebius
end

tol = 1e16;
if any(abs(points)>tol)
  error('Your points are too Herculean in magnitude...please redefine my tolerance');
end

sample_points = [0, 1, -1, 2:100];
N_sample_points = length(sample_points);
q = 1;

done = false;

while q<=N_sample_points
  % Try mapping Inf to samples_points(q)
  H = map_infinity_to(sample_points(q));
  temp = moebius(points, H);
  if all(abs(temp)<tol)
    done = true;
    break
  end
end

if done
  points = temp;
  map = H;
else
  error('I ran out of all the sample points....');
end
