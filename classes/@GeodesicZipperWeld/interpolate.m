function[varargout] = interpolate(self, varargin)
% interpolate -- Interpolates welding map
%
% theta_exterior = interpolate(self, theta_interior)
%
%     Pushes theta_interior values forward through the welding map.
%
% theta_interior = interpolate(self, [], theta_exterior)
%
%     Pulls back theta_exterior_values through the map inverse.
%
% [exterior, interior] = interpolate(self, theta_interior, theta_exterior)
%
%     Simultaneously pushes forward theta_interior ---> exterior,
%                        pulls back interior <--- theta_exterior

%opt = strict_inputs({'theta_int', 'theta_ext'}, {[], []}, [], varargin{:});
if length(varargin)==1
  varargin{2} = [];
end

% Miscellaneous preprocessing
interior_size = size(varargin{1});
varargin{1} = varargin{1}(:); 

exterior_size = size(varargin{2});
varargin{2} = varargin{2}(:);

% First figure out if we're doing interior to exterior or vice versa, etc.
map_to_exterior = not(isempty(varargin{1}));
map_to_interior = not(isempty(varargin{2}));

% If nothing to be done...return
if not(map_to_interior | map_to_exterior)
  varargout{1:nargout} = [];
  return
end
% Append 0 values cuz for some reason I can't figure out how to make this work
% without 0 as a first value.
%varargin{1} = [0; varargin{1}];
%varargin{2} = [0; varargin{2}];
N_interior = length(varargin{1});
N_exterior = length(varargin{2});

% Make exceptions for complex power functions (goes into self.inverse_terminal_map)
slit_interior_limbo = (varargin{1} > self.interior_vertices(end)) & ...
                      (varargin{1} < self.interior_vertices(1));
slit_exterior_limbo = (varargin{1} > self.exterior_vertices(end)) & ...
                      (varargin{1} < self.exterior_vertices(1));

% Invert the Moebius alignment
%[z_interior, z_exterior] = inverse_moebius_alignment(exp(i*opt.theta_int), exp(i*opt.theta_ext), mapdata);
[z_interior, z_exterior] = self.inverse_moebius_alignment(exp(i*varargin{1}), exp(i*varargin{2}));

%[vertices_int, vertices_ext] = inverse_moebius_alignment(mapdata.vertices_in, mapdata.vertices_out, mapdata);
[vertices_int, vertices_ext] = self.inverse_moebius_alignment(... 
         self.interior_disc_vertices, self.exterior_disc_vertices);

% Fix real-valued crap and find which vertices the interpolated points lie
% between
vertices_int = real(vertices_int); vertices_ext = real(vertices_ext);
vertices_int(1) = -Inf; vertices_ext(1) = -Inf;  % By construction of the map

N_vertices = length(vertices_int);

if map_to_exterior
  z_interior = real(z_interior);
  z_interior(abs(mod(varargin{1} - self.interior_interval(1),2*pi))<1e-14) = -Inf;
  [z_interior, sort_order_interior] = sort(z_interior); 

  [bin_count_int, bin_id_int] = histc(z_interior, [vertices_int; Inf]);
  bin_count_int(1) = bin_count_int(1) + bin_count_int(end);
  bin_count_int(end) = [];
  bin_id_int(bin_id_int==(N_vertices+1)) = 1;
else
  bin_id_int = [];
end

if map_to_interior
  z_exterior = real(z_exterior);
  %z_exterior(abs(mod(varargin{2} - self.exterior_interval(1),2*pi))<1e-14) = -Inf;
  [z_exterior, sort_order_exterior] = sort(z_exterior); 

  [bin_count_ext, bin_id_ext] = histc(z_exterior, [vertices_ext; Inf]);
  bin_count_ext(1) = bin_count_ext(1) + bin_count_ext(end);
  bin_count_ext(end) = [];
  bin_id_ext(bin_id_ext==(N_vertices+1)) = 1;
else
  bin_id_ext = [];
end

% Find last (really the first) tooth we must zipup to:
%max_tooth = min([bin_id_int; bin_id_ext]);
max_tooth = min([bin_id_int-1; bin_id_ext-1]);
max_tooth = max([max_tooth 1]);

z = [z_interior; z_exterior];
null_ind = false(size(z));
null_ind = [bin_id_int==N_vertices; bin_id_ext==N_vertices];
interior_ind = null_ind; interior_ind(1:N_interior) = true;
interior_ind(null_ind) = false;
exterior_ind = null_ind; exterior_ind(N_interior+1:end) = true;
exterior_ind(null_ind) = false;
slit_interior_limbo = interior_ind & [slit_interior_limbo; false(size(bin_id_ext))];
slit_exterior_limbo = exterior_ind & [slit_exterior_limbo; false(size(bin_id_ext))];

% Invert the terminal map
z = self.inverse_terminal_map(z, null_ind, interior_ind, exterior_ind, slit_interior_limbo, slit_exterior_limbo);

%slide = select_slider(mapdata.type);
%slider/geodesic_slider

slit_interior = interior_ind;
slit_exterior = exterior_ind;

N_interior_tabled = 0;
N_exterior_tabled = 0;

% I wrote this once upon a clairvoyant time...it still seems to work.
for q = self.N_teeth:-1:max_tooth
  if map_to_exterior
    interior_ind2 = N_interior - N_interior_tabled;
    interior_ind1 = N_interior + 1 - N_interior_tabled - bin_count_int(q+2);
    slit_interior(interior_ind1:interior_ind2) = false;
    N_interior_tabled = N_interior_tabled + bin_count_int(q+2);
  end
  if map_to_interior
    exterior_ind2 = N_interior + N_exterior - N_exterior_tabled;
    exterior_ind1 = N_interior + N_exterior+1 -N_exterior_tabled - bin_count_ext(q+2);
    slit_exterior(exterior_ind1:exterior_ind2) = false;
    N_exterior_tabled = N_exterior_tabled + bin_count_ext(q+2);
  end

  [z] = self.slider('zipup', q, z, null_ind, slit_interior, slit_exterior);
end

% Switch signs for those on bin_count_int/ext(max_tooth)
%switch_bin = 1;
switch_bin = max_tooth;
if map_to_exterior
  i2 = N_interior - N_interior_tabled - bin_count_int(switch_bin+1);;
  i1 = N_interior + 1 - N_interior_tabled - bin_count_int(switch_bin+1) - bin_count_int(switch_bin);
  z(i1:i2) = -real(z(i1:i2));
  null_ind(1:bin_count_int(switch_bin)) = true;  % 'interior' points along z_0 -- z_1
end
if map_to_interior
  i2 = N_interior + N_exterior - N_exterior_tabled - bin_count_ext(switch_bin+1);
  i1 = N_interior + N_exterior+1 - N_exterior_tabled - bin_count_ext(switch_bin+1) - bin_count_ext(switch_bin);
  z(i1:i2) = -real(z(i1:i2));
  null_ind(N_interior+1:N_interior+bin_count_ext(switch_bin)) = true; % 'exterior' points
end

slit_interior(null_ind) = false;
slit_exterior(null_ind) = false;

for q = (max_tooth):self.N_teeth
  [z] = self.slider('unzip', q, z, null_ind, slit_exterior, slit_interior);

  % Need to include new points that were tabled, as well as enforce correct sign
  if map_to_interior
    z(bin_id_ext==(q)) = -abs(real(z(bin_id_ext==(q))));
    N_exterior_tabled = N_exterior_tabled - bin_count_ext(q+2);
    null_ind(slit_exterior) = true;
    slit_exterior(:) = false;
    exterior_ind2 = N_interior + N_exterior - N_exterior_tabled;
    exterior_ind1 = N_interior + N_exterior+1 - N_exterior_tabled - bin_count_ext(q+2);
    slit_exterior(exterior_ind1:exterior_ind2) = true;
  end
  if map_to_exterior
    z(bin_id_int==(q)) = abs(real(z(bin_id_int==(q))));
    N_interior_tabled = N_interior_tabled - bin_count_int(q+2);
    null_ind(slit_interior) = true;
    slit_interior(:) = false;
    interior_ind2 = N_interior - N_interior_tabled;
    interior_ind1 = N_interior + 1 - N_interior_tabled - bin_count_int(q+2);
    slit_interior(interior_ind1:interior_ind2) = true;
  end

end

% Redo the terminal map
z = self.terminal_map(z);

if map_to_interior
  z_exterior = z(N_interior+1:end);
  z_exterior(sort_order_exterior) = z_exterior;
end
if map_to_exterior
  z_interior = z(1:N_interior);
  z_interior(sort_order_interior) = z_interior;
end

% Redo the moebius alignment map
[z_exterior, z_interior] = self.moebius_alignment(z_exterior, z_interior);

% We're done
tempi = self.interior_wrap(angle(z_interior));
tempe = self.exterior_wrap(angle(z_exterior));
% First get rid of first points
%tempi(1) = []; tempe(1) = [];
phi_out = reshape(tempi, interior_size);  % These have been mapped to the exterior
phi_in = reshape(tempe, exterior_size);   % These have been mapped to the interior

% Deal output(s)
if xor(map_to_interior, map_to_exterior)

  if map_to_interior
    varargout{1} = phi_in;
  else
    varargout{1} = phi_out;
  end
  varargout{2} = [];
else
  varargout{1} = phi_out;
  varargout{2} = phi_in;
end
