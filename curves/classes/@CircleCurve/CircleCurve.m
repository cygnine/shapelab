classdef CircleCurve
% CircleCurve -- Object implementation of a circle in 2D space
%
% obj = CircleCurve({a=[], b=[], c=[], z_points, x0=[], r=[]})
%
%     Constructs a circle object. You may either 
%          1.) Specify a, b, and c, three points that the circle passes through,
%          2.) Specify both r and x0. r is the radius of the circle, and x0 is
%              the center
%          3.) Specify a 3-tuple array with a complex entry indicating the
%              points in the plane to interpolate.

properties 
  r = 0;
  x0 = [0; 0];
  m = [];
  isaline = false;
end

methods

  function self = CircleCurve(varargin)
  % CircleCurve -- Class constructor
    persistent strict_inputs interpolate_circle
    if isempty(strict_inputs)
      from labtools import strict_inputs
      from shapelab.curves import interpolate_circle
    end

    inputs = {'a', 'b', 'c', 'r', 'x0','z_points'};
    defaults = {[],[],[],[],[],[]};

    opt = strict_inputs(inputs, defaults, [], varargin{:});

    if isempty(opt.z_points)  % No complex points present
      if isempty(opt.a) | isempty(opt.b)  | isempty(opt.c) % Three points aren't present
        if isempty(opt.r) | isempty(opt.x0)  % Nothing present -- default
          self.r = 0;
          self.x0 = 0;
        else                    % radius + center present
          self.r = opt.r;
          self.x0 = opt.x0;
        end
      else            % three points present
        [self.r, self.x0] = interpolate_circle(opt.a,opt.b,opt.c);
        if isinf(self.r) % Then it's a line
          self.x0 = opt.a;
        end
      end
    else  % 3 complex points present
      realpoints = real(opt.z_points);
      imagpoints = imag(opt.z_points);
      [self.r, self.x0] = interpolate_circle([realpoints(1) imagpoints(1)], ...
                                             [realpoints(2) imagpoints(2)], ...
                                             [realpoints(3) imagpoints(3)]);
      if isinf(self.r) % Then it's a line
        self.x0 = [realpoints(1) imagpoints(1)];
      end
    end
  end

  [] = plot(self);

end

end
