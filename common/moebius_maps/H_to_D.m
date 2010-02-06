function[H] = D_to_H()
% D_to_H -- Maps the disc to the upper-half plane
%
% H = D_to_H()
%
%     The base Moebius map taking the disc D to the upper-half plane H. This
%     base map has the following images:
%
%        0   ------> 1
%        1   ------> +i
%        Inf ------> -1
%
%    If you want to specify points on this mapping, use specify_points instead.

persistent H
if isempty(H)
  H = [-1, i; ...
        1, i];
end
