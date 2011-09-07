function[curves] = init__()
% init__ -- Initializes the curves package of shapelab
%
% nodes = init__()

%curves = recurse_files(pwd);
%pwd_addpath('classes');

curves.module_list = {};
curves.recurse_files = true;
curves.addpaths = {'classes'};
