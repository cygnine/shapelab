function[z] = initial_map(self, z)
% initial_map -- The initial map for the geodesic zipper-type algorithm.
%
% z = initial_map(self, z)
%
%     Performs the initial map for zipper-type algorithms. 

z = i*sqrt(self.moebius_maps.initial_map(z));
