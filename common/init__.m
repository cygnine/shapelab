function[common] = init__()
% init__ -- Initialization file for shapelab/common package
%
% [nodes] = init__()

module_list = {'moebius_maps'};
%common = recurse_files(pwd, module_list);

common.module_list = module_list;
common.recurse_files = true;
common.addpaths = {};
