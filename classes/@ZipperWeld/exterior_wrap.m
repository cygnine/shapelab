function[out] = exterior_wrap(self,inp)
% exterior_wrap -- Wraps real-valued angles to 2*pi-length interval for weld
%
% out = exterior_wrap(self, inp)
%
%     Expresses inp as the values in the interval self.exterior_interval by
%     adding an appropriate multiple of 2*pi.

out = self.exterior_interval(1) + mod(inp - self.exterior_interval(1), 2*pi);
