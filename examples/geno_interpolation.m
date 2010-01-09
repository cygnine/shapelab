% script to test geno interpolation

clear 
close all

from shapelab.geno import geno_interpolant
from shapelab.curves import eval_curvature_samples

load mpeg7_contours;
temp = mpeg7_contour{167};
N = 100;
z = temp(round(linspace(1,length(temp),N+1))); 

z = temp(1:20:(end-6));

[k,s,z0] = geno_interpolant(z, 'closed', true, 'k', 3);

N = 100;

points = zeros([N length(z)]);

for q = 1:length(z);
  s_temp = linspace(s(1,q), s(2,q), N);
  points(:,q) = z0(q) + eval_curvature_samples(k(:,q), s_temp);
end

plot(points, 'r-'); hold on; plot(z, 'bo'); axis equal
