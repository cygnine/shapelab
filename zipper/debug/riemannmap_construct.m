% A script file to test formation of zipper maps.

clear
close all
from shapelab.test_shapes import ellipse
from shapelab.zipper import riemann_maps

Nz = 100;
z = ellipse(Nz);

[int_map, ext_map] = riemann_maps(z);
