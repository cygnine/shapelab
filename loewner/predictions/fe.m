function[out_data] = fe(data)
% fe -- Forward-Euler Loewner evolution prediction
%
% out_data = fe(data)
%
%     Predicts the driving function lambda(s) using the next data point to be
%     `zipped' down in a Loewner evolution. This prediction comes in the form of
%     a value of derivative of lambda and for a stepsize. 

dl_prediction = imag(data.a)/imag(data.g);
  
ds_prediction = real((dl_prediction^2 + 2*data.lambda*dl_prediction - ...
                      2*data.g*dl_prediction + data.lambda^2 + ...
                      data.g^2 - 2*data.lambda*data.g)/(-4));

out_data.ds = ds_prediction;
out_data.dlambda = dl_prediction./ds_prediction;
