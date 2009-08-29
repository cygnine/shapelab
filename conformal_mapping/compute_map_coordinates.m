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
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%          conformal mapping", 2006.

global handles;
inputs = {'z_in', 'w_in', 'z_out', 'w_out', 'winding_number',...
          'zip_magnitude','type'};
defaults = {false, false, false, false, 1, 0.85, 'geodesic'};
opt = handles.common.InputSchema(inputs,defaults,[],varargin{:});
shapelab = handles.shapelab;
moebius = shapelab.common.moebius;
moebius_inv = shapelab.common.moebius_inverse;
switch lower(opt.type)
case 'geodesic'
  fa = shapelab.conformal_mapping.geodesic.base_conformal_map;
case 'slit'
  error('not yet implemented');
case 'zipper'
  error('not yet implemented');
end

z_n = z_n(:);
N = length(z_n);
assert(N>2,'Error: you must give at least three points defining a shape');

%% Append z_out, z_in to end of vector
% If z_out was not specified, assume it's infinity
if isa(opt.z_out,'logical') & opt.z_out==false
  z_n(end+1) = Inf;
else
  z_n(end+1) = opt.z_out;
end

% No matter what z_in is, tack it on
z_n(end+1) = opt.z_in;

%% These are needed as output for the initial map
z_initial = [z_n(1); z_n(2)];

%% After this line, array order is now z(2),z(3),...,z(n),z_out,z_in,z(1)
zeta_n = circshift(z_n,-1);
% This is weird to me: taking the sqrt branch cut at pi instead of 2*pi. I don't
% understand...just following Marshall. If I use a sqrt branch cut at 2*pi (as
% is done in the rest of the algorithm), I shouldn't need the i factor.
m_initial = [1 -z_initial(2);...
             1 -z_initial(1)];
zeta_n = i*sqrt(moebius(zeta_n(2:end), m_initial)); % don't need original z(2)
zeta_n(end) = Inf; % sadly, sqrt(complex Inf) = NaN

unzipped_in = [Inf;0]; % original z_1, z_2
unzipped_out = [Inf;0];

fa_opt.cut_magnitude = opt.zip_magnitude;

for q = 1:(N-2)
  % The next map is defined by the parameter:
  a_array(q) = zeta_n(1);

  %% Apply the map to all the points
  % interior points + infinity + normalization point + original z(1)
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

% Finally, the terminal map. By now, only 3 points are left: z_out, z_in, z(1)
a_array(end+1) = 1/zeta_n(3);
m = [ 1              0; ...
     -a_array(end) 1];
%% Care about where z_out + z_in points go:
zeta_n = -sign(opt.winding_number)*moebius(zeta_n(1:2), m).^2;
% map unzipped points as well:
unzipped_in = -sign(opt.winding_number)*moebius(unzipped_in,m).^2;
unzipped_out = -sign(opt.winding_number)*moebius(unzipped_out,m).^2;

% Some moebius maps:
m1 = [-1, i; ...
       1, i];  % to unit circle
m2 = [0, 1;...
      1 0]; % inversion so unit disc ----> outside unit disc
m4 = [0, -1;...
      -1 0]; % invert m2
m5 = [i, -i; ...
      -1, -1]; % invert m1

% Did the user specify the exterior point?
if isa(opt.w_out,'logical') & opt.z_out==false
  % Map infinity to infinity
  map_to = 0;
else
  map_to = 1/opt.w_out;
end

a = moebius(zeta_n(1), m2*m1); % image of z_out to disc
m3a = [abs(a)/a -abs(a); ...
      -conj(a) 1];  % map a (z_out) to 0
if map_to~=0
  m3b = [1 abs(map_to); ...
    conj(map_to) abs(map_to)/map_to];  % map 0 (z_out) to map_to
else
  m3b = eye(2);
end

m_out = m5*m4*m3b*m3a*m2*m1;
% This thing maps the real line to the real line...therefore normalize it to get
% rid of imaginary (machine eps crap) stuff
m_out = real(m_out/m_out(1,1));

a = moebius(zeta_n(2),m1);
m7 = [abs(a)/a -abs(a); ...
      -conj(a) 1];  % map normalization point to 0
m_in = m5*m7*m1;
m_in = real(m_in/m_in(1,1));

% Total map
if isa(opt.z_in, 'logical') & isa(opt.w_in, 'logical')
  % Don't really do anything
  m_in = eye(2);
else
  unzipped_in = moebius(unzipped_in,m_in);
end

% No matter what, we do the exterior map
unzipped_out = moebius(unzipped_out,m_out);

% Finally, map to unit circle
m = [-1, i;...
      1, i];

unzipped_in = moebius(unzipped_in, m);
unzipped_out = moebius(unzipped_out, m);

% Classes, anyone?
[mapdata.a_array, mapdata.zip_magnitude, mapdata.z_in, mapdata.w_in,...
 mapdata.z_out, mapdata.w_out, mapdata.vertices_in, mapdata.vertices_out,...
 mapdata.winding_number, mapdata.m_in, mapdata.m_out, mapdata.m_initial,...
 mapdata.type] = ...
 deal(a_array, opt.zip_magnitude,opt.z_in, opt.w_in, opt.z_out, opt.w_out, ...
      unzipped_in, unzipped_out, opt.winding_number, m_in, m_out, m_initial,...
      opt.type);
