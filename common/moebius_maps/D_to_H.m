function[H] = D_to_H(varargin)
% D_to_H -- Maps the disc to the upper-half plane
%
% H = D_to_H()
%
%     The base Moebius map taking the disc D to the upper-half plane H. This
%     base map has the following images:
%
%        1  ------> 0
%        +i ------> 1
%        -1 ------> Inf
%
%    If you want to specify points on this mapping, use specify_points instead.

persistent H
if isempty(H)
  H = [i, -i;...
      -1, -1];
end
