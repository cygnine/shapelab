classdef LineCurve
% LineCurve -- Object implementation of a line in 2D space
%
% obj = LineCurve({a=[], b=[], m=[], x0=[]})
%
%     Constructs a line object. You may either 
%          1.) Specify both a and b, two points that the line passes through,
%          2.) Specify both m and x0. m is the slope of the line and x0 is a
%              point that the line passes through.

properties 
  m = 0;
  x0 = [0; 0];
end

methods

  function self = LineCurve(varargin)
  % LineCurve -- Class constructor
    persistent strict_inputs interpolate_line
    if isempty(strict_inputs)
      from labtools import strict_inputs
      from shapelab.curves import interpolate_line
    end

    inputs = {'a', 'b', 'm', 'x0'};
    defaults = {[],[],[],[]};

    opt = strict_inputs(inputs, defaults, [], varargin{:});

    if isempty(opt.a) | isempty(opt.b)  % Two points aren't present
      if isempty(opt.m) | isempty(opt.x0)  % Nothing present -- default
        self.m = 0;
        self.x0 = [0;0];
      else                    % slopt+point present
        self.m = opt.m;
        self.x0 = opt.x0;
      end

    else            % two points present
      [self.m, self.x0] = interpolate_line(opt.a,opt.b);
    end

  end

end

end
