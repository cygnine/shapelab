classdef GeodesicZipperMap < ZipperMap
  properties
  end
  methods
    function self = GeodesicZipperMap(z, varargin)
    % GeodesicZipperMap -- A 'geodesic' zipper-type conformal map
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
