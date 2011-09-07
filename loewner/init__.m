function[loewner] = init__()
% init__ -- Initializes the loewner submodule
%
% [nodes] = init__()

module_list = {'predictions', 'solutions'};
%loewner = recurse_files(pwd, module_list);

loewner.module_list = module_list;
loewner.recurse_files = true;
loewner.addpaths = {};
