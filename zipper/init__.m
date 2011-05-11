function[zipper] = init__()
% init__ -- Initialization file for shapelab/zipper package
%
% [nodes] = init__()

module_list = {'drivers', 'sliders'};

zipper = recurse_files(pwd, module_list);

%zipper.drivers = matlab_import('drivers');
%zipper.sliders = matlab_import('sliders');
