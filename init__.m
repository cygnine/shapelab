function[shapelab] = init__()
% init__ -- Initializes the shapelab package
%
% [nodes] = init__()

shapelab.common = matlab_import('common');
shapelab.conformal_mapping = matlab_import('conformal_mapping');
shapelab.test_shapes = matlab_import('test_shapes');
shapelab.welding = matlab_import('welding');
