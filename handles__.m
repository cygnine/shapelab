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
hs.conformal_mapping.debug.base = fullfile(hs.conformal_mapping.base, 'debug');
hs.conformal_mapping.zipper.base = fullfile(hs.conformal_mapping.base, 'zipper');
  hs.conformal_mapping.zipper.geodesic.base = ...
    fullfile(hs.conformal_mapping.zipper.base, 'geodesic');
  hs.conformal_mapping.zipper.slit.base = ...
    fullfile(hs.conformal_mapping.zipper.base, 'slit');
  hs.conformal_mapping.zipper.zipper.base = ...
    fullfile(hs.conformal_mapping.zipper.base, 'zipper');
hs.test_shapes.base = fullfile(hs.base, 'test_shapes');
hs.welding.base = fullfile(hs.base, 'welding');
hs.examples.base = fullfile(hs.base, 'examples');
hs.common.base = fullfile(hs.base, 'common');

pathadditions = cell(0);
