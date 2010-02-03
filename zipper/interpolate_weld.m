function[phi_out, phi_in] = interpolate_weld(mapdata, varargin)
% interpolate_weld -- Interpolates a welding map
%
% [phi_out, phi_in] = interpolate_weld(mapdata, {theta_in=[], theta_out=[]})
%
%     Interpolates the welding map specified by mapdata. At this stage, only
%     `geodesic'-type maps are supported. The points theta_in are transferred
%     via the weld to phi_out, and the points theta_out are transferred to
%     phi_in.

persistent strict_inputs zip moebius 
persistent inverse_terminal_map assert_enough_points select_slider inverse_moebius_alignment
persistent terminal_map moebius_alignment

if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.common import moebius 
  from shapelab.zipper import inverse_initial_map inverse_terminal_map assert_enough_points
  from shapelab.zipper import select_slider inverse_moebius_alignment
  from shapelab.zipper import terminal_map moebius_alignment
end

opt = strict_inputs({'theta_in', 'theta_out'}, {[], []}, [], varargin{:});

% Miscellaneous preprocessing
map_to_exterior=true; exterior_size = size(opt.theta_out); 
opt.theta_out = opt.theta_out(:); N_exterior = length(opt.theta_out);
map_to_interior=true; interior_size = size(opt.theta_in); 
opt.theta_in = opt.theta_in(:); N_interior = length(opt.theta_in);

% First figure out if we're doing interior to exterior or vice versa, etc.
if isempty(opt.theta_in)
  map_to_exterior=false;
end
if isempty(opt.theta_out)
  map_to_interior=false;
end

% If nothing to be done...return
if not(map_to_interior | map_to_exterior)
  return
end

% Invert the Moebius alignment
[z_interior, z_exterior] = inverse_moebius_alignment(exp(i*opt.theta_in), exp(i*opt.theta_out), mapdata);
[vertices_int, vertices_ext] = inverse_moebius_alignment(mapdata.vertices_in, mapdata.vertices_out, mapdata);


% Fix real-valued crap and find which vertices the interpolated points lie
% between
vertices_int = real(vertices_int); vertices_ext = real(vertices_ext);
vertices_int(1) = -Inf; vertices_ext(1) = -Inf;  % By construction of the map
N_vertices = length(vertices_int);

if map_to_exterior
  z_interior = real(z_interior);
  z_interior(abs(mod(opt.theta_in,2*pi))<1e-14) = -Inf;
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
  z_exterior(abs(mod(opt.theta_out,2*pi))<1e-14) = -Inf;
  [z_exterior, sort_order_exterior] = sort(z_exterior); 
  [bin_count_ext, bin_id_ext] = histc(z_exterior, [vertices_ext; Inf]);
  bin_count_ext(1) = bin_count_ext(1) + bin_count_ext(end);
  bin_count_ext(end) = [];
  bin_id_ext(bin_id_ext==(N_vertices+1)) = 1;
else
  bin_id_ext = [];
end

% Find last (really the first) tooth we must zipup to:
max_tooth = min([bin_id_int; bin_id_ext]);

z = [z_interior; z_exterior];
null_ind = false(size(z));
null_ind(bin_id_ext==N_vertices) = true;
null_ind(bin_id_int==N_vertices) = true;
interior_ind = null_ind; interior_ind(1:N_interior) = true;
interior_ind(null_ind) = false;
exterior_ind = null_ind; exterior_ind(N_interior+1:end) = true;
exterior_ind(null_ind) = false;

% Invert the terminal map
z = inverse_terminal_map(z, mapdata, null_ind, interior_ind, exterior_ind);

slide = select_slider(mapdata.type);

slit_interior = interior_ind;
slit_exterior = exterior_ind;

N_interior_tabled = 0;
N_exterior_tabled = 0;

for q = mapdata.N_teeth:-1:max_tooth
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

  [z] = slide('zipup', q, z, mapdata, null_ind, slit_interior, slit_exterior);
end

% Switch signs for those on bin_count_int/ext(1)
if map_to_exterior
  i2 = N_interior - N_interior_tabled - bin_count_int(2);;
  i1 = N_interior + 1 - N_interior_tabled - bin_count_int(2) - bin_count_int(1);
  z(i1:i2) = -real(z(i1:i2));
  null_ind(1:bin_count_int(1)) = true;  % 'interior' points along z_0 -- z_1
end
if map_to_interior
  i2 = N_interior + N_exterior - N_exterior_tabled - bin_count_ext(2);
  i1 = N_interior + N_exterior+1 - N_exterior_tabled - bin_count_ext(2) - bin_count_ext(1);
  z(i1:i2) = -real(z(i1:i2));
  null_ind(N_interior+1:N_interior+bin_count_ext(1)) = true; % 'exterior' points
end

slit_interior(null_ind) = false;
slit_exterior(null_ind) = false;

for q = (max_tooth):mapdata.N_teeth
  [z] = slide('unzip', q, z, mapdata, null_ind, slit_exterior, slit_interior);
  if map_to_interior
    N_exterior_tabled = N_exterior_tabled - bin_count_ext(q+2);
    null_ind(slit_exterior) = true;
    slit_exterior(:) = false;
    exterior_ind2 = N_interior + N_exterior - N_exterior_tabled;
    exterior_ind1 = N_interior + N_exterior+1 - N_exterior_tabled - bin_count_ext(q+2);
    slit_exterior(exterior_ind1:exterior_ind2) = true;
  end
  if map_to_exterior
    N_interior_tabled = N_interior_tabled - bin_count_int(q+2);
    null_ind(slit_interior) = true;
    slit_interior(:) = false;
    interior_ind2 = N_interior - N_interior_tabled;
    interior_ind1 = N_interior + 1 - N_interior_tabled - bin_count_int(q+2);
    slit_interior(interior_ind1:interior_ind2) = true;
  end

end

% Redo the terminal map
z = terminal_map(z, mapdata);

if map_to_interior
  z_exterior = z(N_interior+1:end);
  z_exterior(sort_order_exterior) = z_exterior;
end
if map_to_exterior
  z_interior = z(1:N_interior);
  z_interior(sort_order_interior) = z_interior;
end

% Redo the moebius alignment map
[z_exterior, z_interior] = moebius_alignment(z_exterior, z_interior, mapdata);

% We're done
phi_out = angle(z_interior);  % These have been mapped to the exterior
phi_in = angle(z_exterior);   % These have been mapped to the interior
