% script to test geno interpolation

clear 
close all



from shapelab.geno import geno_interpolant
from shapelab.curves import curvature_coefficients as curv_coeffs
from shapelab.curves import extend_curve
from shapelab.curves import circle_coefficients eval_curvature_samples

cd ../../test_shapes
load mpeg7_contours;
cd ../geno/debug;
%temp = mpeg7_contour{1355};
temp = mpeg7_contour{167};
N = 100;
z = temp(round(linspace(1,length(temp),N+1))); 

z = temp(1:20:(end-6));

%theta = linspace(0, 2*pi, 30);
%theta(end) = [];
%z = exp(i*theta);
%z = z + 0.1*randn(size(z)) + 0.1*i*rand(size(z));

[k,s,z0] = geno_interpolant(z, 'closed', true, 'k', 3);

N = 100;

points = zeros([N length(z)]);

for q = 1:length(z);
  s_temp = linspace(s(1,q), s(2,q), N);
  points(:,q) = z0(q) + eval_curvature_samples(k(:,q), s_temp);
end
