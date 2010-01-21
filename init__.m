function[shapelab] = init__()
% init__ -- Initializes the shapelab package
%
% [nodes] = init__()

shapelab.common = matlab_import('common');
shapelab.zipper = matlab_import('zipper');
shapelab.test_shapes = matlab_import('test_shapes');
shapelab.welding = matlab_import('welding');
shapelab.loewner = matlab_import('loewner');
shapelab.curves = matlab_import('curves');
shapelab.geno = matlab_import('geno');
shapelab.extensions = matlab_import('extensions');
shapelab.wp = matlab_import('wp');
