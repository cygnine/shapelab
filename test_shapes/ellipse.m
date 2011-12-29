function[z] = ellipse(N, varargin)
% ellipse -- Sample points from an ellipse
%
% z = ellipse(N, {a=1, b=0.75, rotation=0, centroid=0})
%
%     Takes N complex-valued samples from an ellipse centered at (0,0) with
%     major axis length 2*a and minor axis length 2*b. The major axis is located
%     along the x-axis -- this can be changed by giving a radian rotation. The
%     center location can be augmented by giving a complex scalar for
%     'centroid'. 
%
%     The N samples are parametrically equispaced, where the parameterization is
%     given by
%
%       z = centroid + exp(i*rotation)*(a*cos(t) + i*b*sin(t)),
%
%     for t in [0, 2*pi].

persistent input_parser parser
if isempty(input_parser)
  from labtools import input_parser

  inputs = {'a', 'b', 'rotation', 'centroid'};
  defaults = {1, 0.75, 0, 0};
  [opt,parser] = input_parser(inputs, defaults, [], varargin{:});
else
  parser.parse(varargin{:});
  opt = parser.Results;
end

t = linspace(0, 2*pi, N+1).';
t(end) = [];

z = opt.centroid + exp(i*opt.rotation)*(opt.a*cos(t) + i*opt.b*sin(t));
