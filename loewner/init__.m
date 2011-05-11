function[loewner] = init__()
% init__ -- Initializes the loewner submodule
%
% [nodes] = init__()

module_list = {'predictions', 'solutions'};

loewner = recurse_files(pwd, module_list);
%loewner.predictions = matlab_import('predictions');
%loewner.solutions = matlab_import('solutions');
