function[shapelab] = init__()
% init__ -- Initializes the shapelab package
%
% [nodes] = init__()

module_list = {'common', 'zipper', 'test_shapes', 'welding', ...
               'loewner', 'curves', 'geno', 'extensions', ...
               'wp', 'diffs1', 'teichmuller'};

shapelab = recurse_files(pwd, module_list);

pwd_addpath('classes');

%shapelab.common = matlab_import('common');
%shapelab.zipper = matlab_import('zipper');
%shapelab.test_shapes = matlab_import('test_shapes');
%shapelab.welding = matlab_import('welding');
%shapelab.loewner = matlab_import('loewner');
%shapelab.curves = matlab_import('curves');
%shapelab.geno = matlab_import('geno');
%shapelab.extensions = matlab_import('extensions');
%shapelab.wp = matlab_import('wp');
%shapelab.diffs1 = matlab_import('diffs1');
%shapelab.teichmuller = matlab_import('teichmuller');
