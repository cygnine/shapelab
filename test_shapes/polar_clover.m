function[z] = polar_clover(N, varargin)
% [z] = polar_clover(N, {lobes=4, lobe_depth_ratio=0.5, lobe_max_radius=1, 
%                        node_locations='boundary', M=25})
%
%     Constructs a nodal representation of a polar-coordinate clover with
%     'lobes' lobes. lobe_depth_ratio represents the 'intensity' of the lobes.
%     For value 0, there is no lobe; when set to 1 all the interiors of the
%     lobes collapse to the origin. node_locations can take on values
%     'boundary', 'interior', and 'exterior'. For 'interior' and 'exterior', the
%     optional value M specifies the resolution in the R direction as given in
%     shapelab/common/polar_linspace. For 'exterior', points up to 2 x
%     lobe_max_radius are used. lobe_max_radius determines the maximal value of
%     |z| this function returns. 

global packages;
inputs = {'lobes', 'lobe_depth_ratio', 'node_locations', 'M', 'lobe_max_radius'};
defaults = {4, 0.5, 'boundary', 25, 1};
opt = packages.labtools.input_schema(inputs, defaults, [], varargin{:});

if strcmpi(opt.node_locations, 'boundary')
  theta = linspace(0,2*pi,N+1);
  theta = theta(1:N).';
  r = ones(size(theta));
elseif strcmpi(opt.node_locations, 'interior')
  linopt.r0 = 0;
  linopt.r1 = (opt.M-1)/opt.M;
  z = packages.shapelab.common.polar_linspace(opt.M,N,linopt);
  theta = angle(z);
  r = abs(z);
elseif strcmpi(opt.node_locations, 'exterior')
  linopt.r1 = 2;
  linopt.r0 = 1+1/opt.M;
  z = packages.shapelab.common.polar_linspace(opt.M,N,linopt);
  theta = angle(z);
  r = abs(z);
end

r_min = opt.lobe_max_radius*(1-opt.lobe_depth_ratio);
r_lobe = (opt.lobe_max_radius*opt.lobe_depth_ratio)/2;

z = r.*((r_min + r_lobe) + r_lobe*cos(theta*opt.lobes)).*exp(i*theta);
