function[zipper] = init__()
% init__ -- Initialization file for shapelab/zipper package
%
% [nodes] = init__()

module_list = {'drivers', 'sliders'};
%zipper = recurse_files(pwd, module_list);

zipper.module_list = module_list;
zipper.recurse_files = true;
zipper.addpaths = {};
