function[mapdata] = compute_map_coordinates(z_n,varargin)
% [z_initial,a_array,zeta_n] = compute_map_coordinates(z_n,{z_in=false,z_out=false,winding_number=1})
%
%     The main workhorse routine for computing a conformal map from the points
%     z_n to points on the unit circle. 
%
%     Output mapdata contains:
%            a_array
%            a_cut_bias
%            z_initial
%            zeta_n
%            normalization
%
%     Computes an array of complex values a_array, with length n-2 where n is
%     the length of the input vector z_n. Also returns z_initial=[z_n(1),
%     z_n(2)], and zeta_n, the location of the final node. These data
%     corresponds to the parameters defining the composition maps of the
%     `geodesic' algorithm for conformal mapping. See [1].
%
%     The initial mapping sending z_n(1) off to infinity, and the terminal
%     mapping sending the real line to the unit circle are assumed to be
%     understood without reference. See implementation of evaluate_map.
%
%     The optional arguments z_in and z_out normalize the map so that f(z_in) =
%     z_out. If you try to force the map to do an unnatural thing like mapping a
%     point inside the convex hull to somewhere outside the unit disc, it will
%     do this without complaining. If you give z_in without giving z_out or vice
%     versa, it ignores normalization. winding_number determines the orientation
%     of the input data with respect to the shape.
%
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%          conformal mapping", 2006.

global handles;
inputs = {'z_in', 'z_out', 'winding_number'};
defaults = {[],[],1};
opt = handles.common.InputSchema(inputs,defaults,[],varargin{:});
moebius = handles.shapelab.common.moebius;
fa = handles.shapelab.conformal_mapping.geodesic.base_conformal_map;

z_n = z_n(:);
N = length(z_n);
assert(N>2,'Error: you must give at least three points defining a shape');

a_array = zeros([N-3,1]);
a_cut_bias = true([N-3,1]);
normalization_persistence = zeros([N-3 1]);

if not(isempty(opt.z_in) && isempty(opt.z_out))
  z_n(end+1) = opt.z_in;
  normalize = true;
else
  normalize = false;
end

% These are needed as output for the initial map
z_initial = [z_n(1); z_n(2)];

% Temporary array, first send z_n(1) to infinity
zeta_n = circshift(z_n,-1);
zeta_n = i*sqrt(moebius(zeta_n(2:end), [1, -z_n(2); 1 -z_n(1)]));
zeta_n(end) = Inf;

%p1 = plot(zeta_n, 'b.'); hold on; p2 = plot(zeta_n(end-1), 'r.');

for q = 1:(N-2)
  % Determine next parameter a
  a_array(q) = zeta_n(1);

  zeta_n = fa(zeta_n(2:end),a_array(q));

  if normalize
    %if q<(N-2);
      %temp = zeta_n(end-1);
      %normalization_persistence(q) = abs(temp)^2/real(temp);
      normalization_persistence(q) = N/max(abs(zeta_n));
      %zeta_n = moebius(zeta_n,[1 0; -1/normalization_persistence(q) 1]);
      zeta_n = zeta_n*normalization_persistence(q);
      %set(p1, 'xdata', real(zeta_n), 'ydata', imag(zeta_n));
      %set(p2, 'xdata', real(zeta_n(end-1)), 'ydata', imag(zeta_n(end-1)));
      %drawnow;
      %pause;
    %end
  end
  zeta_n(end-1)
end

if normalize
  temp = zeta_n;
  zeta_n = zeta_n(2);
  temp = -(moebius(temp, [1 0; -1/zeta_n 1])).^2;

  % Map to unit circle:
  temp = moebius(temp, [-1 i; 1 i]);
  normalization = [temp(1), opt.z_out];
else
  normalization = [0,0];
end

[mapdata.z_initial, mapdata.a_array, mapdata.a_cut_bias, ...
  mapdata.normalization_persistence, mapdata.zeta_n, mapdata.normalization] = ...
  deal(z_initial, a_array, a_cut_bias, normalization_persistence, zeta_n, normalization);
