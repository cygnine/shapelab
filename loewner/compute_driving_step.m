function[new_lambda, dlambda, ds, a, an, g, gn] = compute_driving_step(a, an, g, gn, lambda, q,M);

persistent invert_a
if isempty(invert_a)
  from shapelab.loewner import invert_a_unzip as invert_a
end

pflags = isinf(g);
mflags = isinf(gn);

end_ind = 1 + q;

dlambda = zeros([M 1]);
ds = zeros([M 1]);
new_lambda = zeros([M+1 1]);
new_lambda(1) = lambda;

for qq = 1:M

  dl_prediction = imag(a(end_ind))/imag(g(end_ind));
  
  ds_prediction = real((dl_prediction^2 + 2*lambda*dl_prediction - ...
                        2*g(end_ind)*dl_prediction + lambda^2 + ...
                        g(end_ind)^2 - 2*lambda*g(end_ind))/(-4));

  dlambda(qq) = dl_prediction/ds_prediction;

  % But now only take a step as big as ds_prediction*qq/M
  ds(qq) = ds_prediction*1/(M+1-qq);

  % Evolve lambda:
  new_lambda(qq+1) = new_lambda(qq) + ds(qq)*dlambda(qq);

  % Evolve the rest of the points:
  a = a + ds(qq)*(2 - g.*dlambda(qq));


  % Not yet unzipped
  an(end_ind:end) = a(end_ind:end);
  
  % Negative real axis:
  an(1:(end_ind-1)) = an(1:(end_ind-1)) + ds(qq)*(2 - gn(1:(end_ind-1)).*dlambda(qq));

  lambda = new_lambda(qq+1);

  [g,gn] = invert_a(a, an, lambda);
end

a(end_ind) = -1/2*(new_lambda(end)).^2;
an(end_ind) = -1/2*(new_lambda(end)).^2;
g(end_ind) = new_lambda(end);
gn(end_ind) = new_lambda(end);


%dl = imag(a(q+1))/imag(g(q+1));
%ds = real((dl^2 + 2*lambda*dl - 2*g(q+1)*dl + lambda^2 + g(q+1)^2 - ...
%           2*lambda*g(q+1))/(-4));

% An approximation to derivative of lambda:
%dlambda = dl/ds;

% Evolve lambda:
%lambda = lambda + dl;

% Evolve the rest of the points:
%a = a + ds*(2-g.*dlambda);

%a(q+1) = -1/2*(lambda.^2);  % Explicitly enforce -- machine roundoff crap

% The positive real axis:
%an(q+1:end) = a(q+1:end);

% The negative real axis:
%an(1:q) = an(1:q) + ds*(2-gn(1:q).*dlambda);

% Recover g:
%[g,gn] = invert_a(a, an, lambda);

a(pflags) = Inf;
g(pflags) = Inf;
an(mflags) = Inf;
gn(mflags) = Inf;
