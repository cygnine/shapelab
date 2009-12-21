% This script shows an example for usage of the delta angle denoiser.

clear 
close all

from shapelab.curves import delta_angle_denoise as denoise

load mpeg7_contours
z0 = mpeg7_contour{167}(1:20:end-6);
winterp = 100;

% Let's do four denoisings:
w{1} = denoise(z0, 'k', 3, 'scale', 'dyadic','interpolation_weight',winterp);
  titles{1} = 'Dyadic scaling, $k=3$ fine scales';
w{2} = denoise(z0, 'k', 5, 'scale', 'dyadic','interpolation_weight',winterp);
  titles{2} = 'Dyadic scaling, $k=5$ fine scales';
w{3} = denoise(z0, 'k', 3, 'scale', 'linear','interpolation_weight',winterp);
  titles{3} = 'Linear scaling, $k=3$ fine scales';
w{4} = denoise(z0, 'k', 5, 'scale', 'linear','interpolation_weight',winterp);
  titles{4} = 'Linear scaling, $k=5$ fine scales';

figure;
for q = 1:4;
  subplot(2,2,q); hold on; 
  plot(z0, 'k.');
  plot(z0, 'r--');
  plot(w{q}, 'b');
  set(title(titles{q}), 'interpreter', 'latex');
  axis equal
  axis off
end
