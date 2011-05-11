function[wp] = init__()
% init__ -- Initializes the shapelab.wp package
%
% [nodes] = init__()

module_list = {'teichon'};

wp = recurse_files(pwd, module_list);
%wp.teichon = matlab_import('teichon');
