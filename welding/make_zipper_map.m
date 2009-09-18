function[mapdata] = make_zipper_map(z,varargin)
% [mapdata] = make_zipper_map(z,{type='zipper', shape_0 = false, 
%                                shape_infinity=Inf, 
%                                winding_number=1})
%
%     A wrapper for the zipper-type conformal mapping constructor. z are the
%     vertices of the shape. shape_0 is the location inside that shape that is
%     mapped to 0 inside the unit disc. shape_infinity is the location outside
%     the shape that is mapped to infinity. Only the sign of winding_number is
%     used in the computations. 

global handles;
inputs = {'type', 'shape_0', 'shape_infinity', 'winding_number','visualize'};
defaults = {'zipper', false, Inf, 1,0};
opt = handles.common.input_schema(inputs, defaults, [], varargin{:});

mopt.type = opt.type;
mopt.winding_number = opt.winding_number;
if not(isa(opt.shape_0, 'logical'))
  mopt.z_in = opt.shape_0;
  mopt.w_in = 0;
end

mopt.z_out = opt.shape_infinity;
mopt.w_out = Inf;
mopt.visualize = opt.visualize;

mapdata = handles.shapelab.conformal_mapping.zipper.compute_map_coordinates(z, mopt);
