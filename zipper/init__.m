function[zipper] = init__()
% init__ -- Initialization file for shapelab/zipper package
%
% [nodes] = init__()

zipper = recurse_files(pwd);
zipper.geodesic = matlab_import('geodesic');
zipper.slit = matlab_import('slit');
zipper.zipper = matlab_import('zipper');
zipper.drivers = matlab_import('drivers');
zipper.sliders = matlab_import('sliders');
