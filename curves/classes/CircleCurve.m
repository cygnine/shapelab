classdef CircleCurve
% CircleCurve -- Object implementation of a circle in 2D space
%
% obj = CircleCurve({a=[], b=[], c=[], x0=[], r=[]})
%
%     Constructs a circle object. You may either 
%          1.) Specify a, b, and c, three points that the circle passes through,
%          2.) Specify both r and x0. r is the radius of the circle, and x0 is
%              the center

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

    inputs = {'a', 'b', 'm', 'x0'};
    defaults = {[],[],[],[]};

    opt = strict_inputs(inputs, defaults, [], varargin{:});

    if isempty(opt.a) | isempty(opt.b)  | isempty(opt.c) % Three points aren't present
      if isempty(opt.r) | isempty(opt.x0)  % Nothing present -- default
        self.r = 0;
        self.x0 = [0;0];
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

  end

end

end
