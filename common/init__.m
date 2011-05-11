function[common] = init__()
% init__ -- Initialization file for shapelab/common package
%
% [nodes] = init__()

module_list = {'moebius_maps'};

common = recurse_files(pwd, module_list);

%common.moebius_maps = matlab_import('moebius_maps');
