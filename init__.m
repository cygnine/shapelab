function[shapelab] = init__()
% init__ -- Initializes the shapelab package
%
% [nodes] = init__()

module_list = {'common', 'zipper', 'test_shapes', 'welding', ...
               'loewner', 'curves', 'geno', 'extensions', ...
               'wp', 'diffs1', 'teichmuller'};
%shapelab = recurse_files(pwd, module_list);
%pwd_addpath('classes');

shapelab.module_list = module_list;
shapelab.recurse_files = true;
shapelab.addpaths = {'classes'};
