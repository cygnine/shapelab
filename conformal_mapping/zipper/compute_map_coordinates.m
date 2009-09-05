function[mapdata] = compute_map_coordinates(z_n,varargin)
% [z_initial,a_array,zeta_n] = compute_map_coordinates(z_n,{z_in=false,w_in=false,
%                 z_out=false, w_out=false, winding_number=1,
%                 zip_magnitude=0.85, type='geodesic'})
%
%     The main workhorse routine for computing a conformal map from the points
%     z_n to points on the unit circle. 
%
%     Output mapdata contains:
%            a_array
%            zip_magnitude
%            z_in
%            w_in
%            z_out
%            w_out
%            vertices_in
%            vertices_out
%            type
%
%     The output mapdata contains all the necessary information for travelling
%     back and forth on the conformal map defined by the shape's vertices z_n
%     and their corresponding inner and outer values on the unit circle,
%     vertices_in and vertices_out.
%
%     This is the main work file for zipper-type algorithms from [1]. The 'base'
%     map may be interchanged to form the 'geodesic', 'slit', and 'zipper'
%     algorithms.
%    
%     The map takes z_in (inside the original shape) to w_in (inside the unit
%     circle). If z_in is specified without w_in, then w_in is assumed to be 0.
%     If neither z_in nor w_in is specified, no interior Moebius transform is
%     performed. 
%
%     The map takes z_out (outside the original shape) to w_out (outside the
%     unit circle.) If z_out is specified without w_out, then w_out is assumed
%     to be infinity. If neither z_out nor w_out is specified, the map takes
%     infinity to infinity.
%
%     The winding_number determines the orientation of the input data with
%     respect to the shape.
%
%     The option zip_magnitude refers to how each 'base' map is normalized.
%     Empirically, values O(1) produce stable results.
%
%     The `type' of mapping can be: 'geodesic', 'slit', or 'zipper'.
%
%     A note on normalization: this function normalizes the interior and
%     exterior moebius maps so that the first vertex of the shape lies at the
%     mapped point z=1 (both interior and exterior).
%
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%          conformal mapping", 2006.

global handles;
inputs = {'z_in', 'w_in', 'z_out', 'w_out', 'winding_number',...
          'zip_magnitude','type'};
defaults = {false, false, Inf, Inf, 1, 0.85, 'geodesic'};
opt = handles.common.InputSchema(inputs,defaults,[],varargin{:});
shapelab = handles.shapelab;
zip = shapelab.conformal_mapping.zipper;
moebius = shapelab.common.moebius;
csqrt = handles.shapelab.common.positive_angle_square_root;
moebius_inv = shapelab.common.moebius_inverse;
switch lower(opt.type)
case 'geodesic'
  fa = zip.geodesic.base_conformal_map;
  zipper = false;
case 'slit'
  fa = zip.slit.base_conformal_map;
  zipper = false;
case 'zipper'
  fa = zip.zipper.base_conformal_map;
  zipper = true;
case 'zipper_weld'
  fa = zip.zipper.base_conformal_map;
  fa_geo = zip.geodesic.base_conformal_map;
  zipper = true;
otherwise
  error(['Unrecognized algorithm specification ''' opt.type '''']);
end

% Some initial computations, making sure we have enough points
z_n = z_n(:);
N = length(z_n);
assert(N>2,'Error: you must give at least three points defining a shape');
if strcmpi(opt.type, 'zipper')
  assert(mod(N,2)==0, 'Error: the zipper algorithm requires an even number of points');
end

%% Append z_out, z_in to end of vector
% If z_out was not specified, assume it's infinity
z_n(end+1) = opt.z_out;

% No matter what z_in is, tack it on
z_n(end+1) = opt.z_in;

%% These are needed as output for the initial map
if zipper
  z_initial = [z_n(1); z_n(2); z_n(3)];
else
  z_initial = [z_n(1); z_n(2)];
end

%% Initial map:
% After this line, array order is now z(2),z(3),...,z(n),z_out,z_in,z(1)
zeta_n = circshift(z_n,-1);

if zipper
  m_initial = [(z_initial(2) - z_initial(1))*[1, -z_initial(3)]; ...
               (z_initial(2) - z_initial(3))*[1, -z_initial(1)]];
  zeta_n = csqrt(moebius(zeta_n(3:end), m_initial));
  zeta_n(end) = Inf; % sadly, sqrt(complex Inf) = NaN

  unzipped_in = [Inf; -1; 0]; % original z_1, z_2, z_3
  unzipped_out = [Inf; 1; 0]; % original z_1, z_2, z_3
else
  m_initial = [1 -z_initial(2);...
               1 -z_initial(1)];
  zeta_n = i*csqrt(moebius(zeta_n(2:end), m_initial)); % don't need original z(2)
  zeta_n(end) = Inf; % sadly, sqrt(complex Inf) = NaN

  unzipped_in = [Inf;0]; % original z_1, z_2
  unzipped_out = [Inf;0];
end
%%

%% Initialization for looping over teeth
fa_opt.cut_magnitude = opt.zip_magnitude;
if zipper
  N_teeth = (N-2)/2 - 1;
else
  N_teeth = N-2;
end

a_array = zeros([N_teeth+1,1]);
c_array = zeros([N_teeth+1,1]); % Only needed for zipper, really
%%

%% Looping over teeth
for q = 1:N_teeth
  if zipper
    % Need next two points to define the map:
    c_array(q) = zeta_n(1);
    a_array(q) = zeta_n(2);

    % Apply the map to all the points
    % interior points + infinity + z_in + original z(1)
    temp = zeros(size(zeta_n(3:end)));
    temp(end) = 2;
    fa_opt.point_id = temp;
    zeta_n = fa(zeta_n(3:end),c_array(q),a_array(q), fa_opt);

    % append new points to unzipped_in, unzipped_out
    unzipped_in(end+1:end+2) = [c_array(q); a_array(q)];
    unzipped_out(end+1:end+2) = [c_array(q); a_array(q)];
    temp = 2*ones(size(unzipped_in), 'int8'); temp(end-2:end) = 1;
    fa_opt.point_id = temp;
    [unzipped_in,garbage] = fa(unzipped_in,c_array(q),a_array(q),fa_opt);
    [garbage,unzipped_out] = fa(unzipped_out,c_array(q),a_array(q),fa_opt);

    % we know:
    unzipped_in(end) = 0; unzipped_out(end) = 0;
  else
    % The next map is defined by the parameter:
    a_array(q) = zeta_n(1);

    % Apply the map to all the points
    % interior points + infinity + z_in + original z(1)
    temp = zeros(size(zeta_n(2:end)));
    temp(end) = 2;
    fa_opt.point_id = temp;
    zeta_n = fa(zeta_n(2:end),a_array(q),fa_opt);

    temp = 2*ones(size(unzipped_in), 'int8'); temp(end) = 1;
    fa_opt.point_id = temp;
    [unzipped_in,garbage] = fa(unzipped_in,a_array(q),fa_opt);
    [garbage,unzipped_out] = fa(unzipped_out,a_array(q),fa_opt);

    % Now add 0 to the list of zipped_in/out points:
    unzipped_in(end+1) = 0;
    unzipped_out(end+1) = 0;
  end
end
%%

%% Terminal map
if strcmpi(opt.type, 'zipper')
  % Finally, the terminal map. By now, only 4 points are left: z(N), z_out, z_in, z(1)
  c_array(end) = zeta_n(1);
  a_array(end) = 1./zeta_n(4);
  m = [1 0; ...
       -a_array(end) 1];
  alpha = pi - angle(moebius(c_array(end),m));
  zeta_n = sign(opt.winding_number)*(exp(-i*(pi-alpha))*moebius(zeta_n(2:3),m)).^(pi/alpha);

  % map unzipped points as well:
  unzipped_in = sign(opt.winding_number)*(exp(-i*(pi-alpha))*moebius(unzipped_in,m)).^(pi/alpha);
  unzipped_out = sign(opt.winding_number)*(exp(-i*(pi-alpha))*moebius(unzipped_out,m)).^(pi/alpha);
elseif strcmpi(opt.type, 'zipper_weld')
  % Use the last point z(N) to do a geodesic arc, then finish off in same way as
  % geodesic algorithm
  c_array(end) = zeta_n(1);

  temp = zeros(size(zeta_n(2:end)));
  temp(end) = 2;
  fa_opt.point_id = temp;
  zeta_n = fa_geo(zeta_n(2:end),c_array(end),fa_opt);

  temp = 2*ones(size(unzipped_in), 'int8'); temp(end) = 1;
  fa_opt.point_id = temp;
  [unzipped_in,garbage] = fa_geo(unzipped_in,c_array(end),fa_opt);
  [garbage,unzipped_out] = fa_geo(unzipped_out,c_array(end),fa_opt);

  % Now add 0 to the list of zipped_in/out points:
  unzipped_in(end+1) = 0;
  unzipped_out(end+1) = 0;

  % Now do same terminal map as in geodesic algorithm
  a_array(end) = 1./zeta_n(3);
  m = [ 1              0; ...
       -a_array(end) 1];
  %% Care about where z_out + z_in points go:
  zeta_n = -sign(opt.winding_number)*moebius(zeta_n(1:2), m).^2;

  % map unzipped points as well:
  unzipped_in = -sign(opt.winding_number)*moebius(unzipped_in,m).^2;
  unzipped_out = -sign(opt.winding_number)*moebius(unzipped_out,m).^2;
else
  % Finally, the terminal map. By now, only 3 points are left: z_out, z_in, z(1)
  a_array(end) = 1/zeta_n(3);
  m = [ 1              0; ...
       -a_array(end) 1];
  %% Care about where z_out + z_in points go:
  zeta_n = -sign(opt.winding_number)*moebius(zeta_n(1:2), m).^2;

  % map unzipped points as well:
  unzipped_in = -sign(opt.winding_number)*moebius(unzipped_in,m).^2;
  unzipped_out = -sign(opt.winding_number)*moebius(unzipped_out,m).^2;
end
%%

%% PSL(2)-type maps
% Some moebius maps:
m1 = [-1, i; ...
       1, i];  % to unit circle
m2 = [0, 1;...
      1 0]; % inversion so unit disc ----> outside unit disc
m4 = [0, -1;...
      -1 0]; % invert m2
m5 = [i, -i; ...
      -1, -1]; % invert m1

% Always specify the exterior map
map_to = 1/opt.w_out;
%%if isa(opt.w_out,'logical') & opt.z_out==false
%if isa(opt.w_out,'logical')
%  % Map infinity to infinity
%  opt.w_out = Inf;
%  map_to = 0;
%else
%  map_to = 1/opt.w_out;
%end

a = moebius(zeta_n(1), m2*m1); % image of z_out to disc
m3a = [abs(a)/a -abs(a); ...
      -conj(a) 1];  % map a (z_out) to 0
if map_to~=0
  m3b = [1 abs(map_to); ...
    conj(map_to) abs(map_to)/map_to];  % map 0 (z_out) to map_to
else
  m3b = eye(2);
end

z_1 = moebius(unzipped_out(1), m3b*m3a*m2*m1);
m_rotate_out = [exp(-i*angle(z_1)), 0; 0, 1]; % Just rotate this point to z=1

m_out = m5*m4*m_rotate_out*m3b*m3a*m2*m1;
% This thing maps the real line to the real line...therefore normalize it to get
% rid of imaginary (machine eps crap) stuff
m_out = real(m_out/max(m_out(:)));

a = moebius(zeta_n(2),m1);
m7 = [abs(a)/a -abs(a); ...
      -conj(a) 1];  % map normalization point to 0
z_1 = moebius(unzipped_in(1), m7*m1);
m_rotate_in = [exp(-i*angle(z_1)), 0; 0, 1];
m_in = m5*m_rotate_in*m7*m1;
m_in = real(m_in/max(m_in(:)));

% Total map
if isa(opt.z_in, 'logical') & isa(opt.w_in, 'logical')
  % Don't really do anything
  m_in = eye(2);
else
  unzipped_in = moebius(unzipped_in,m_in);
end

% No matter what, we do the exterior map
unzipped_out = moebius(unzipped_out,m_out);
%%

%% End up on the unit circle
% Finally, map to unit circle
m = [-1, i;...
      1, i];

unzipped_in = moebius(unzipped_in, m);
unzipped_out = moebius(unzipped_out, m);

% Hmmm...screaming new data structure
[mapdata.a_array, mapdata.zip_magnitude, mapdata.z_in, mapdata.w_in,...
 mapdata.z_out, mapdata.w_out, mapdata.vertices_in, mapdata.vertices_out,...
 mapdata.winding_number, mapdata.m_in, mapdata.m_out, mapdata.m_initial,...
 mapdata.type, mapdata.c_array, mapdata.N_teeth] = ...
 deal(a_array, opt.zip_magnitude,opt.z_in, opt.w_in, opt.z_out, opt.w_out, ...
      unzipped_in, unzipped_out, opt.winding_number, m_in, m_out, m_initial,...
      opt.type, c_array, N_teeth);
