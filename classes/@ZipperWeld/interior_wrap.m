function[out] = interior_wrap(self,inp)
% interior_wrap -- Wraps real-valued angles to 2*pi-length interval for weld
%
% out = interior_wrap(self, inp)
%
%     Expresses inp as the values in the interval self.interior_interval by
%     adding an appropriate multiple of 2*pi.

out = self.interior_interval(1) + mod(inp - self.interior_interval(1), 2*pi);
