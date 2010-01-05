function[mapdata, otherdata, z] = loewner_driver(z, z_in, z_out, Q, visdata)
% loewner_driver -- The main driver for the Loewner evolution unzipping
%
% [mapdata, otherdata, z] = loewner_driver(z, z_in, z_out, Q, visdata)

persistent visualize loewner_step compute_a
if isempty(visualize)
  from shapelab.zipper.drivers import visualize
  from shapelab.loewner import compute_driving_step as loewner_step
  from shapelab.loewner import compute_a
end

%z_in = [z_in; zeros([Q 1])];
%z_out = [z_out; zeros([Q 1])];

M = 30; % Number of sub-samples

lambda = zeros([Q*M+1 1]);
dlambda = zeros([Q*M 1]);

N = length(z_in);

g = [z_out; z];
gn = [z_in; z];

lambda(1) = g(N);
a = compute_a(g, lambda(1));
an = compute_a(gn, lambda(1));
s = zeros([Q*M+1 1]);

for q = 1:Q

    ind1 = 2 + (q-1)*M;
    ind2 = 1 + q*M;

    [lambdan, dlambdan, ds, a, an, g, gn] = loewner_step(a, an, g, gn, lambda(ind1-1), q+N-1, M);

    s(ind1:ind2) = s(ind1-1) + cumsum(ds);
    lambda(ind1:ind2) = lambdan(2:end);
    dlambda((ind1-1):(ind2-1)) = dlambdan;

    if visdata.visualize
      visualize(g((q+N):(end-3)), g(1:(q+1)), gn(1:(q+1)), visdata)
    end
end

z = g(end-2:end);

[otherdata.lambda, otherdata.dlambda, otherdata.s] = deal(...
           lambda,           dlambda,           s);

[mapdata.unzipped_out, mapdata.unzipped_in] = deal(...
         g(1:(end-3)),        gn(1:(end-3)));
