function[zipper] = init__()
% init__ -- Initialization file for shapelab/zipper package
%
% [nodes] = init__()

zipper = recurse_files(pwd);
zipper.drivers = matlab_import('drivers');
zipper.sliders = matlab_import('sliders');
