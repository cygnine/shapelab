function[z]= evaluate_inverse_map(w,mapdata,varargin)
% [z]= evaluate_inverse_map(w,mapdata,{point_id=zeros(size(w))})
%
%     Input mapdata contains:
%            a_array
%            z_initial
%            zeta_n
%            normalization
%
%     Evaluates the inverse of the zipper-type conformal map defined by mapdata.
%     The optional input 'point_id' should be specified as follows:
%
%     z(point_id==0) -----> (default) points somewhere inside or outside the
%     unit circle that you want to map inside or outside the shape.
%
%     z(point_id==1) -----> points on the interior boundary of the unit disc
%     that you want to map onto the shape. 
%
%     z(point_id==2) -----> points on the exterior boundary of the unit disc
%     that you want to map onto the shape.
%
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%     conformal mapping", 2006.

global handles;
opt = handles.common.InputSchema({'point_id'}, {zeros(size(w))}, [], varargin{:});
ifa = handles.shapelab.conformal_mapping.geodesic.inverse_base_conformal_map;
moebius = handles.shapelab.common.moebius;
moebius_inv = handles.shapelab.common.moebius_inverse;
csqrt = handles.shapelab.common.positive_angle_square_root;

%[z_initial, a_array, a_cut_bias, normalization_persistence, zeta_n, normalization] = ...
%  deal(mapdata.z_initial, mapdata.a_array, mapdata.a_cut_bias, ...
%    mapdata.normalization_persistence, mapdata.zeta_n, mapdata.normalization);

% We must deal with three cases:
w0 = w(opt.point_id==0);
w1 = w(opt.point_id==1);
w2 = w(opt.point_id==2);

z = zeros(size(w));

%% Hokay, here we go

% If we're on the interior/exterior proper, no hocus-pocus is necessary
if any(opt.point_id==0)
  % Meh, just in case
  if any(abs(w0)==1)
    error('Point(s) specified as strictly not on unit circle looks like it is');
  end
  ext_points = abs(w0)>1;
  int_points = abs(w0)<1;

  % map unit circle to R:
  m = [-1, i;...
        1, i];
  w0 = moebius_inv(w0,m);

  % apply inv(m_in) on interior
  w0(int_points) = moebius_inv(w0(int_points), mapdata.m_in);

  % apply inv(m_out) on exterior
  w0(ext_points) = moebius_inv(w0(ext_points), mapdata.m_out);

  % This should map H (interior) ----> left-half plane (right-half if winding
  % number is negative)
  w0 = csqrt(-sign(mapdata.winding_number)*w0);

  % Invert the terminal map (moebius)
  m = [1, 0;...
       -mapdata.a_array(end), 1];
  w0 = moebius_inv(w0,m);

  % Now loop through a_array
  N = length(mapdata.a_array)+1;

  ifa_opt.temp = zeros(size(w0));
  ifa_opt.cut_magnitude = mapdata.zip_magnitude;
  for q = (N-2):-1:1
    w0 = ifa(w0, mapdata.a_array(q), ifa_opt);
  end

  % Now it's easy-peasy:
  w0 = moebius_inv((w0/i).^2, mapdata.m_initial);
  z(opt.point_id==0) = w0;
end

if any(opt.point_id==1)
  error('Not coded yet');
end

if any(opt.point_id==2)
  error('Not coded yet');
end

%% Un-normalize
%w = dab(z, normalization(2), normalization(1));
%
%% Map from unit circle to half-plane
%w = moebius(w, [i -i; -1 -1]);
%if opt.boundary
%  w = real(w);
%end
%
%% Invert terminal map
%w = csqrt(mapdata.terminal_sign*w,'cut_bias', false);
%w = moebius(w, [1 0; 1/zeta_n 1]);
%
%% Invert sequential maps
%for q = length(a_array):-1:1
%  %if q<=length(a_array) && normalization_persistence(q)~=0
%  %  w = w/normalization_persistence(q);
%  %end
%  w = ifa(w,a_array(q));
%end
%
%% Invert initial map
%w = (w/i).^2;
%w = moebius(w, [-z_initial(1) z_initial(2); -1 1]);
