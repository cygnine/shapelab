function[out_data] = lserk4(data)
% lserk4 -- Low-storage explicit RK-4 Loewner evolution prediction
%
% out_data = lserk4(data)
%
%     Predicts the driving function lambda(s) using the next data point to be
%     `zipped' down in a Loewner evolution. This prediction comes in the form of
%     a value of derivative of lambda and for a stepsize. 

persistent lserk_temp rk4 fe
if isempty(rk4)
  from odesolve.coeffs import lserk4 as lserk_temp
  rk4 = lserk_temp();
  from shapelab.loewner.predictions import fe
end

dl_prediction = imag(data.a)/imag(data.g);
  
ds_prediction = real((dl_prediction^2 + 2*data.lambda*dl_prediction - ...
                      2*data.g*dl_prediction + data.lambda^2 + ...
                      data.g^2 - 2*data.lambda*data.g)/(-4));

out_data.ds = ds_prediction;
out_data.dlambda = dl_prediction./ds_prediction;
