classdef ZipperMap < ConformalMap
  properties 
    z_in 
    w_in
    type
    N
    z
    a_array
    w
    moebius_maps
  end

  methods
    function self = ZipperMap(z, varargin)
    % ZipperMap -- A zipper-type RiemannMap.
    %
    % obj = ZipperMap(z, 
    %     This class defines a type of RiemannMap w = f(z) that is a discrete
    %     solution to the Loewner evolution equation. The method of construction is
    %     sequential iteration over the values z. The input z is an ordered
    %     collection of complex-value points that lie on the boundary of some
    %     connected region U in the complex plane. The output is then a map that
    %     taking U to some other region (usually the unit disc).
    %
    %     The constructor of this class is not meant to be an end-user tool --
    %     instances are generated from shapelab.zipper functions.

    persistent strict_inputs
    if isempty(strict_inputs)
      from labtools import strict_inputs
    end

    self = self@RiemannMap()

    inputs = {'type', 'z_in', 'w_in', 'N', 'z', 'w', 'a_array', 'moebius_maps'};
    defaults = {'geodesic', zeros(0), zeros(0), 0, zeros(0), zeros(0), zeros(0), struct()};
    opt = strict_inputs(inputs, defaults, {}, varargin{:});
    fnames = fieldnames(opt);
    for q = 1:length(fnames);
      self = setfield(self, fnames{q}, getfield(opt, fnames{q}));
    end
    self.z = z;

    % Perform some rudimentary sanity checks:
    assert(length(self.z) == length(self.w), ...
           'Point and images vectors z and w must have the same size');
    self.N = length(self.z);

    valid_types = {'geodesic', 'slit', 'zipper', 'loewner'};
    assert(any(strcmpi(self.type, valid_types)), ...
           'You must give a valid zipper ''type''.');
    end

  end

end
