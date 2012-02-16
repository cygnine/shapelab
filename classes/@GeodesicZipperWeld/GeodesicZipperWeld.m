classdef GeodesicZipperWeld < ZipperWeld
% GeodesicZipperlWeld -- A 'geodesic' type of zipper welding map.
%
% self = GeodesicZipperWeld()
%     A functional object representing a 'geodesic' flavor of conformal
%     zipper weld. See [1] for details.
%
%     [1]: D. Marshall, S. Rohde; "Convergence of a variant of the zipper
%     algorithm for conformal mapping", SIAM J. Numerical Analysis, v45, no 6.

  properties(SetAccess=protected)
    a_array
    N
  end
  properties(Access=protected)
    derivative_at_inf; % variable used only temporarily for construction of map
  end

  methods
    function self = GeodesicZipperWeld(z, varargin)

      persistent input_parser parser
      if isempty(parser)
        from labtools import input_parser

        [opt,parser] = input_parser({'tooth_length', 'visualize'}, {0.05, false}, []);
      end

      if isempty(z)
        % User asking us to take discrete weld evaluations
        [newargs] = varargin(3:end);
      else
        % User asking us to take samples on shape
        newargs = varargin;
      end

      parser.parse(newargs{:});
      opt = parser.Results;

      self = self@ZipperWeld(newargs{:});
      self.tooth_length = opt.tooth_length;
      if isempty(z)
        self = self.weld_from_samples(varargin{1}, varargin{2});
        return;
      end

      fd_at_inf = [1e6+2 1e6+3].';

      self.z = z;

      if isempty(self.center_location)
        [N, self.N_teeth, z] = self.assert_enough_points(z, fd_at_inf, self.inf_location, mean(z), z(1));
      else
        [N, self.N_teeth, z] = self.assert_enough_points(z, fd_at_inf, self.inf_location, self.center_location, z(1));
      end
      self.N = N;

      [z, self] = self.calculate_initial_map(z);

      w = z;  % w are 'interior' points, z are 'exterior' points

      visdata = self.initialize_visualization(opt.visualize, z);
      visdata.visualize = opt.visualize;

      %self.a_array = zeros([length(z) - 4 1]);
      %self.moebius_maps.tooth_maps = cell([length(z) - 4 1]);
      self.a_array = zeros([N-1 1]);
      self.moebius_maps.tooth_maps = cell(size(self.a_array));

      for q = 1:self.N_teeth
        [w, z, self] = self.unzip(q, w, z);
        %self.visualize(z(q+2:end-3), w(1:q+2), z(1:q+2), visdata);
        self.visualize(z(q+2:N), w(1:q+2), z(1:q+2), visdata);
      end

      [w, z, self] = self.calculate_terminal_map(w, z);
      [w, z, self] = self.calculate_moebius_alignment(w, z);
      % Done
    end

    [self] = rechart(self, varargin);
    varargout = interpolate(varargin);
    z = map_to_shape(self, z, varargin);
    z = map_from_shape(self, z, varargin);
    %[z_int,z_ext] = separate_to_shape(self, theta_int, theta_ext);
    val = calculate_derivative_at_inf(self,zinf);
    [self] = match(self, other)
  end

  methods(Access=protected)
    [z, self] = calculate_initial_map(self, z)
    z = initial_map(self, z);
    [z_int, z_ext, self] = calculate_inverse_initial_map(self, z_int, z_ext)
    visdata = initialize_visualization(self, vistf, z)
    [w, z, self] = unzip(self, tooth, z, w);
    [z_int, z_ext, self] = zipup(self, tooth, z_int, z_ext);
    z = slider(self, direction, tooth, z, interior, slit_interior, slit_exterior);
    [] = visualize(self, z, z_in, z_out, visdata);
    [w, z, self] = calculate_terminal_map(self, w, z);
    [z_int, z_ext, self] = calculate_inverse_terminal_map(self, z_int, z_ext);
    z = terminal_map(self, z);
    [w, z, self] = calculate_moebius_alignment(self, w, z);
    [z] = inverse_terminal_map(self, z, interior, slit_interior, slit_exterior, slit_interior_limbo, slit_exterior_limbo);
    [zi, ze] = inverse_moebius_alignment(self, zi, ze);
    [z_int, z_ext, self] = calculate_inverse_moebius_alignment(self, z_int, z_ext);
    [N, N_teeth, z] = assert_enough_points(self, z, varargin)
    self = weld_from_samples(self, theta_int, theta_ext);
  end

end
