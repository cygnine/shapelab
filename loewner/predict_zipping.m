function[data] = predict_zipping(data, type)
% predict_zipping -- Computations for prediction of Loewner evolution
%
% data = predict_zipping(data, type)
%
%     Predicts lambda', the derivative of the driving function for Loewner
%     evolutions of slits.

persistent rk2
if isempty(rk2)
  from shapelab.loewner.predictions import rk2
end

switch type
case 'rk2'
  data = rk2(data)
end
