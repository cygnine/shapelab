function[psi_out, psi_in] = interpolate_fingerprint(vertices_int, vertices_ext, varargin)
% interpolate_fingerprint -- Interpolates a fingerprint
%
% [psi_out, psi_in] = interpolate_fingerprint(vertices_int, vertices_ext,
% {theta_int=[], theta_ext=[]})
%
%     Interpolates a fingerprint using the given (vertices_int,
%     vertices_ext) pairs. These are assumed to come from a function over a
%     2*pi-length domain *that is monotonic*. 
%
%     The fingerprint psi satisfies psi(theta+2*pi) = psi(theta) + 2*pi. The
%     fingerprint is defined as 
%
%        psi = inv(phi_ext) o phi_int,
%
%     where phi_int (phi_ext) is a map of the interior (exterior) of the disc to
%     the interior (exterior) of a shape. 
%
%     The input theta_int are angles on the interior of the unit disc and get
%     mapped via the weld to the external values psi_out. theta_ext and psi_in
%     form a similar pair.
%
% Assumptions: 
%
%    - The shape has a positive winding number
%    - No chart computations
%    - Geodesic tooth-length of 0.1

persistent wrap strict_inputs find_moebius moebius moebius_inverse
persistent force_inverse_terminal_map zipup unzip terminal_map fd
if isempty(strict_inputs)
  from labtools import strict_inputs
  from labtools import interval_wrap as wrap
  from shapelab.common.moebius_maps import specify_points as find_moebius
  from shapelab.common import moebius moebius_inverse
  from shapelab.zipper.drivers import force_inverse_terminal_map terminal_map
  from shapelab.loewner.solutions import normal_linear_slit_zip as zipup
  from shapelab.loewner.solutions import normal_linear_slit_unzip as unzip
  from finite_difference import difference_derivative_periodic as fd
end

N_vertices = length(vertices_int);

assumptions.type = 'geodesic';
assumptions.tooth_length = 0.1;
assumptions.winding_number = 1;
assumptions.moebius_maps.interior_terminal = eye(2);
assumptions.moebius_maps.exterior_terminal = eye(2);
assumptions.moebius_map.terminal_map = eye(2);
assumptions.N_teeth = N_vertices - 2;

opt = strict_inputs({'theta_int', 'theta_ext'}, {[], []}, [], varargin{:});
map_to_interior = not(isempty(opt.theta_ext));
N_exterior = length(opt.theta_ext);
map_to_exterior = not(isempty(opt.theta_int));
N_interior = length(opt.theta_int);

vertices_int = vertices_int(:);
vertices_ext = vertices_ext(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Preprocessing: make the first map the `best-resolved' one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find approximations to derivative
dapprox = fd(vertices_int, vertices_ext, 2, [vertices_int(1) ,2*pi + vertices_int(1)]);
% Make point 1 the location with the smallest abs(log(|f'|))
[garbage, ground_zero] = min(abs(log(abs(dapprox))));

% Hokay, rotate everything to ground zero
vertices_int = circshift(vertices_int, 1-ground_zero);
vertices_ext = circshift(vertices_ext, 1-ground_zero);

interval = [0, 2*pi];
fprint_int_bias = vertices_int(1);
fprint_ext_bias = vertices_ext(1);
vertices_int = wrap(vertices_int, interval);
vertices_ext = wrap(vertices_ext, interval);

psi_out_size = size(opt.theta_int);
psi_in_size = size(opt.theta_ext);
opt.theta_int = wrap(opt.theta_int(:), interval);
opt.theta_ext = wrap(opt.theta_ext(:), interval);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Moebius alignment here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do alignment for input into terminal zipper map:
% Send: vertices(end)   ---> 0
%       vertices(1)     ---> Inf
%       vertices(end-1) ---> somewhere on real line

% Find point somewhere on other side:
middle_index = mod(ground_zero + floor(N_vertices/2) - 1, N_vertices) + 1;

% Perform alignment: wrap fingerprint so that vertices(1) ---> (0,0)
vertices_int = exp(i*vertices_int);
vertices_ext = exp(i*vertices_ext);
opt.theta_int = exp(i*opt.theta_int);
opt.theta_ext = exp(i*opt.theta_ext);
H_int = eye(2);
H_ext = eye(2);

middle_point_int = vertices_int(middle_index);
middle_point_ext = vertices_ext(middle_index);

H_int = find_moebius([vertices_int(1) middle_point_int, vertices_int(end)], ...
                     [Inf,            -1,                0]);

H_ext = find_moebius([vertices_ext(1) vertices_ext(end) middle_point_ext], ...
                     [Inf,            0,                -1]);


vertices_int = real(moebius(vertices_int, H_int));
opt.theta_int = real(moebius(opt.theta_int, H_int));
vertices_ext = real(moebius(vertices_ext, H_ext));
opt.theta_ext = real(moebius(opt.theta_ext, H_ext));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get ready for zippering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vertices_int(1) = -Inf;
vertices_ext(1) = -Inf;
vertices_int(end) = 0;
vertices_ext(end) = 0;

tol = 1e-14;
opt.theta_int(abs(opt.theta_int)<tol) = Inf;
opt.theta_ext(abs(opt.theta_ext)<tol) = Inf;
% Do you want to make sure z(supposed_to_be_inf) = -Inf?
[opt.theta_int, sort_order_interior] = sort(opt.theta_int);
[opt.theta_ext, sort_order_exterior] = sort(opt.theta_ext);

if map_to_exterior
  [bin_count_int, bin_id_int] = histc(opt.theta_int, [vertices_int; Inf]);
  bin_count_int(1) = bin_count_int(1) + bin_count_int(end);
  bin_count_int(end) = [];
  bin_id_int(bin_id_int==(N_vertices+1))=1;
else
  bin_id_int = [];
end
if map_to_interior
  [bin_count_ext, bin_id_ext] = histc(opt.theta_ext, [vertices_ext; Inf]);
  bin_count_ext(1) = bin_count_ext(1) + bin_count_ext(end);
  bin_count_ext(end) = [];
  bin_id_ext(bin_id_ext==(N_vertices+1)) = 1;
else
  bin_id_ext = [];
end

% bin_id == q -----> point lies on slit between vertices(q) and vertices(q+1)
% For a geodesic map, there are (N_vertices - 2) slits -- the first two vertices
% are connected via a circular arc, as are the last vertex and the first vertex.

z = [opt.theta_int; opt.theta_ext];

% Now form matrices that have indexing information about various bins
if map_to_exterior
  bin_indices_int = zeros([N_vertices 2]);
  bin_indices_int(:,2) = cumsum(bin_count_int(:));
  bin_indices_int(1,1) = 1;
  bin_indices_int(2:end,1) = bin_indices_int(1:end-1,2)+1;
  pick_bin_int = logical(spalloc(length(z), N_vertices, length(z))); 
  for q = 1:N_vertices
    is = bin_indices_int(q,1):bin_indices_int(q,2);
    pick_bin_int(is,q) = true;
  end
else
  bin_indices_int = [];
  pick_bin_int = logical(spalloc(length(z), N_vertices, length(z)));
end
if map_to_interior
  bin_indices_ext = zeros([N_vertices 2]);
  bin_indices_ext(:,2) = cumsum(bin_count_ext(:));
  bin_indices_ext(1,1) = 1;
  bin_indices_ext(2:end,1) = bin_indices_ext(1:end-1,2)+1;
  bin_indices_ext = bin_indices_ext + N_interior;
  pick_bin_ext = logical(spalloc(length(z), N_vertices, length(z))); 
  for q = 1:N_vertices
    is = bin_indices_ext(q,1):bin_indices_ext(q,2);
    pick_bin_ext(is,q) = true;
  end
else
  bin_indices_ext = [];
  pick_bin_ext = logical(spalloc(length(z), N_vertices, length(z)));
end

% Now, e.g. z(bin_indices_int(34,1):bin_indices_int(34,2)) are all the points
% that are in bin 34 on the interior. So are z(pick_bin_int(34,:)).

% Find last (really the first) tooth we must zipup to:
max_tooth = min([bin_id_int; bin_id_ext]);


slit_interior = [true([N_interior 1]); false([N_exterior 1])];
slit_exterior = not(slit_interior);

interior = pick_bin_int(:,end) | pick_bin_ext(:,end);
slit_interior = slit_interior & not(interior);
slit_exterior = slit_exterior & not(interior);

% For the terminal map:
% interior == anything in interior of the upper half-plane (i.e. on the
%             vertices(1) -- vertices(end) slit
% slit_interior == anything on the negative real axis
% slit_exterior == anything on the positive real axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now push into terminal zipper map:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[vertices_int, vertices_ext, z, assumptions] = ...
  force_inverse_terminal_map(vertices_int, vertices_ext, z,...
  assumptions, interior, slit_interior, slit_exterior);

% It turns out it's easier to code this backward zipping if we treat `interior'
% points as if they're on the slit -- since the zipup action doesn't care about
% whether things are on slits or not, this is a benign preference.
interior = false(size(z));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now go through all the slits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho = assumptions.tooth_length;  % meh, alias
M = zeros([assumptions.N_teeth - max_tooth + 1, 2]);
h = eye(2);
for q = assumptions.N_teeth:-1:1
  bin_id = q + 1;

  % First task: compute the Moebius map that 'symmetrizes' the next vertices.
  z1 = vertices_ext(q+1);
  z2 = vertices_int(q+1);
  M(q,:) = (inv([z1*rho, rho; -z2*rho, -rho])*[z1; z2]).'; % Faster than calling specify_points
  % The full map:
  h(2,:) = M(q,:);

  % Now do the following:
  %  1.) apply the previous Moebius map to all affected points
  %  2.) Take away points that are about to be unzipped -- make them their negative
  if map_to_exterior
    z(slit_interior) = moebius(z(slit_interior), h);

    is = pick_bin_int(:,bin_id);
    slit_interior(is) = false;
    z(is) = -z(is);
  end
  if map_to_interior
    z(slit_exterior) = moebius(z(slit_exterior), h);
    is = pick_bin_ext(:,bin_id);
    slit_exterior(is) = false;
    z(is) = -z(is);
  end

  is = slit_interior | slit_exterior;
  % Now apply the zipup map to the remaining points
  z(is) = zipup(i*rho, z(is));

  % Of course, no matter what we need to apply the moebius map and zipup the
  % remaining vertices
  vertices_int(1:q+1) = zipup(i*rho, moebius(vertices_int(1:q+1), h));
  vertices_ext(1:q+1) = zipup(i*rho, moebius(vertices_ext(1:q+1), h));
end

% The last bin must be 'symmetrized' and 'negativized'. To do this, map the
% circular arc orthogonal to 0 and vertices(1) to the upper-half of imaginary axis.
h1 = find_moebius([0 vertices_int(1)/2+i*vertices_int(1)/2 vertices_int(1)], [0 i Inf]); % vertices_int(1) and vertices_ext(1) should be identical
% h1 should be real-valued
h1 = real(h1/h1(1,1));
is = pick_bin_int(:,1) | pick_bin_ext(:,1);
z(is) = -moebius(z(is),h1);

% Now begin collecting points for unzipping
interior = is;  % the only interior points now are those in bin 1
z(interior) = moebius_inverse(z(interior), h1);

% Now unzip each slit
for q = 1:assumptions.N_teeth
  bin_id = q;

  % Things should be symmetric: unzip
  z(interior) = unzip(i*rho, z(interior));

  % Godd%*@! these multivalued functions
  z(pick_bin_int(:,q)) = abs(z(pick_bin_int(:,q)));
  z(pick_bin_ext(:,q)) = -abs(z(pick_bin_ext(:,q)));

  % Add new bin points to interior:
  interior = interior | pick_bin_int(:,q+1) | pick_bin_ext(:,q+1);

  % Undo the Moebius map
  h(2,:) = M(q,:);
  z(interior) = moebius_inverse(z(interior), h);
end

% Apply terminal map to everyone
z = real(terminal_map(z, assumptions));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi_out = wrap(angle(moebius_inverse(z(1:N_interior), H_ext)), interval);
psi_in = wrap(angle(moebius_inverse(z(N_interior+1:end), H_int)), interval);

psi_out(sort_order_interior) = psi_out;
psi_in(sort_order_exterior) = psi_in;

psi_out = reshape(psi_out, psi_out_size);
psi_in = reshape(psi_in, psi_in_size);
