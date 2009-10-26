function[conformal_mapping] = init__()
% init__ -- Initialization file for shapelab/conformal_mapping package
%
% [nodes] = init__()

conformal_mapping = recurse_files(pwd);
conformal_mapping.zipper = matlab_import('zipper');
