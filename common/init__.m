function[common] = init__()
% init__ -- Initialization file for shapelab/common package
%
% [nodes] = init__()

common = recurse_files(pwd);
common.moebius_maps = matlab_import('moebius_maps');
