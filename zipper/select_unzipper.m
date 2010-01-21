function[fcn] = select_unzipper(type)
% select_unzipper -- Selects zipper/unzipper for zipper-type algorithms
%
% fcn = select_unzipper(type)

switch lower(type)
case 'geodesic'
  from shapelab.zipper.sliders import geodesic_unzipper as fcn
case 'slit'
  from shapelab.zipper.sliders import slit_unzipper as fcn
case 'loewner'
  from shapelab.zipper.sliders import loewner_unzipper as fcn
otherwise
  error('Unrecognized zipper type');
end
