function[fcn] = select_slider(type)
% select_slider -- Selects zipper/unzipper for zipper-type algorithms
%
% fcn = select_slider(type)

switch lower(type)
case 'geodesic'
  from shapelab.zipper.sliders import geodesic_slider as fcn
case 'slit'
  from shapelab.zipper.sliders import slit_slider as fcn
case 'loewner'
  from shapelab.zipper.sliders import loewner_slider as fcn
otherwise
  error('Unrecognized zipper type');
end
