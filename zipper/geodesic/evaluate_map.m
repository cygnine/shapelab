function[w]= evaluate_map(z,mapdata,varargin)
% [w]= evaluate_map(z,mapdata)
%
%     Input mapdata contains:
%            a_array
%            a_cut_bias
%            z_initial
%            zeta_n
%            normalization
%            terminal_sign
%
%     Evaluates the basic zipper-type conformal map at the points z. If z
%     corresponds points *on* the shape, be forewarned that this program will do
%     some funny things unless the input you give is *exactly* on the shape. 
%     algorithm', specified by the data z_initial, a_array, and zeta_n. See
%     compute_map_coordinates and [1] for details of what these inputs are. The
%     output w is the map evaluated at the input locations z. normalization
%     specifies the normalization to be taken for the conformal map. 
%
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%     conformal mapping", 2006.

persistent input_schema fa moebius dab
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import moebius
  from shapelab.zipper.geodesic import base_conformal_map as fa
  from shapelab.common import disc_a_to_b as dab
end

[z_initial, a_array, a_cut_bias, normalization_persistence, zeta_n, normalization] = ...
  deal(mapdata.z_initial, mapdata.a_array, mapdata.a_cut_bias, ...
    mapdata.normalization_persistence, mapdata.zeta_n, mapdata.normalization);

% Initial map
w = i*sqrt(moebius(z,[1 -z_initial(2); 1 -z_initial(1)]));
w(not(isfinite(w))) = Inf;
%figure; p1 = plot(w, 'b.');

for q = 1:length(a_array)
  w = fa(w,a_array(q));
  %if q<=length(a_array) && normalization_persistence(q)~=0
    %w = moebius(w, [1 0; -1/normalization_persistence(q) 1]);
    %w = w*normalization_persistence(q);
    %set(p1, 'xdata', real(w), 'ydata', imag(w));
    %drawnow;
  %end
end

% Terminal map:
w = mapdata.terminal_sign*(moebius(w,[1 0; -1/zeta_n 1])).^2;

% Map to unit circle:
w = moebius(w, [-1 i; 1 i]);

% Normalize
w = dab(w,normalization(1), normalization(2));
