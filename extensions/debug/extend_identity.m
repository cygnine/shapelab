% A script to make sure that the function written does extend the identity

clear
close all

from shapelab.extensions import *

% Test strictly the upper half-plane
is = 0.1:0.1:10;
x = linspace(-10, 10, 100);

z = zeros([length(x) length(is)]);
for q = 1:length(is);
  z(:,q) = x + i*is(q);
end

h = @(x) sign(x).*sqrt(abs(x));
zext = zeros(size(z));
for q = 1:length(is);
  zext(:,q) = ba_driver(h, z(:,q));
end

plot3(real(z), imag(z), abs(zext), 'r.');
hold on;
plot3(real(z), imag(z), 0*real(z), 'k');
plot3(real(z.'), imag(z.'), 0*real(z.'), 'k');

figure;
mesh(abs(zext));

plot(z, 'r.');
figure;
plot(zext, 'b.');
