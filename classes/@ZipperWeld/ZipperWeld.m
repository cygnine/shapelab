classdef ZipperWeld < ConformalWeld
% ZipperWeld -- A conformal welding map of 'zipper' type.
%
% self = ZipperWeld()
%     A functional object representing a 'zipper' conformal welding map.
%     Zipper-type conformal welds are described in [1].
%
%     [1]: D. Marshall, S. Rohde; "Convergence of a variant of the zipper
%     algorithm for conformal mapping", SIAM J. Numerical Analysis, v45, no 6.

  properties(SetAccess=protected)
    z
    inf_location
    inf_image
    center_location
    center_image
    winding_number
    tooth_length
    M
    moebius_maps
    N_teeth
    interior_vertices
    interior_disc_vertices
    exterior_vertices
    exterior_disc_vertices
    interior_interval = [0, 2*pi];
    exterior_interval = [0, 2*pi];
  end
  properties(Access=protected)
    theta_tol = 1e-12;
  end

  methods
    function self = ZipperWeld(varargin)

      persistent input_parser parser
      if isempty(parser)
        from labtools import input_parser

        inputs = {'inf_location', 'inf_image', 'center_location', 'center_image', 'winding_number'};
        defaults = {Inf, Inf, 0, 0, 1};
        [opt,parser] = input_parser(inputs, defaults, [], varargin{:});
      else
        parser.parse(varargin{:});
        opt = parser.Results;
      end

      self.inf_location = opt.inf_location;
      self.inf_image = opt.inf_image;
      self.center_location = opt.center_location;
      self.center_image = opt.center_image;
      self.winding_number = opt.winding_number;

      self.moebius_maps.H_to_D = MoebiusMap([-1, i;...
                                              1, i]);
      self.moebius_maps.invert_D = MoebiusMap([0, 1;...
                                               1, 0]);
      self.moebius_maps.D_to_H = MoebiusMap([i, -i; ...
                                             -1, -1]);
    end

    function [] = plot(self,varargin)
    % plot -- Plots the welding map vertices on [0,2 pi]
    %
    %     Plots the x-y coordinates for interior-exterior vertex coordinates on
    %     S^1.

      if length(varargin) == 0
        plot(sort(self.interior_vertices), sort(self.exterior_vertices), 'b.-');
      else
        plot(sort(self.interior_vertices), sort(self.exterior_vertices), varargin{:});
      end
      axis([self.interior_interval, self.exterior_interval]);
    end


    function[self] = set.interior_interval(self,newint)
      % Effectively ignore newint(2:end)
      self.interior_interval = real([newint(1), newint(1)+2*pi]);
    end
    function[self] = set.exterior_interval(self,newint);
      % Effectively ignore newint(2:end)
      self.exterior_interval = real([newint(1), newint(1)+2*pi]);
    end
    %function[self] = set.interior_vertices(newvals)
    function[self] = set.interior_vertices(self,newvals)
      % Just call wrap
      self.interior_vertices = self.interior_wrap(newvals);
      % This is a fix for a common problem:
      if self.interior_vertices(1) > self.interior_vertices(2);
        self.interior_vertices(1) = self.interior_vertices(1) - 2*pi;
      end

      %self.interior_interval = self.interior_vertices;
    end
    function[self] = set.exterior_vertices(self,newvals)
      % Just call wrap
      self.exterior_vertices = self.exterior_wrap(newvals);
      % This is a fix for a common problem:
      if self.exterior_vertices(1) > self.exterior_vertices(2);
        self.exterior_vertices(1) = self.exterior_vertices(1) - 2*pi;
      end

      %self.exterior_interval = self.exterior_vertices;
    end

    % Externally-defined methods_
    out = interior_wrap(self, inp);
    out = exterior_wrap(self, inp);
  end

  methods(Abstract,Access=protected)
    [N, N_teeth, z] = assert_enough_points(z, varargin)
  end

end
