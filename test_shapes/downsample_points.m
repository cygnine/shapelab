function[z] = downsample_points(z, N)
% downsample_points -- Downsamples a shape
%
% z = downsample_points(z, N)
%
%     Downsamples a collection of points z to N samples. Note that in the
%     languate of signal processing, this does NOT do downsampling. This
%     procedure samples N points from the vector z -- it does not do
%     interpolation. 

z = z(:);
tol = 1e-12;
if abs(z(1)-z(end))<tol
  z(end) = [];
end

M = length(z);
z = [z; z(1)];
if M<N
  fprintf('Warning: You requested more samples than I was given. Doing nothing\n');
  return
end

z = z(round(linspace(1,M+1,N+1)));
z(end) = [];
