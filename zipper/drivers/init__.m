function[drivers] = init__()
% init__ -- Initialization file for shapelab/zipper/drivers package
%
% [drivers] = init__()

drivers = recurse_files(pwd);
