function[hs,pathadditions] = handles__()
% [hs,pathadditions] = handles__()
%
%     Returns directory pointers for common module in hs. pathadditions is a
%     cell array with a string in each element indicated paths to add to the
%     global path structure. 

% This is by default
hs.base = fileparts(mfilename('fullpath'));

% Add subdirectories manually
hs.conformal_mapping.base = fullfile(hs.base, 'conformal_mapping');
  hs.conformal_mapping.geodesic_algorithm.base = ...
    fullfile(hs.conformal_mapping.base, 'geodesic_algorithm');
hs.test_shapes.base = fullfile(hs.base, 'test_shapes');
hs.common.base = fullfile(hs.base, 'common');

pathadditions = cell(0);
