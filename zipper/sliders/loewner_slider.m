function[w_n, z_n, mapdata] = loewner_slider(direction, tooth, z_n, w_n, mapdata)
% loewner_slider -- Unzips one tooth using a predictive Loewner slit evolution
%
% [w_n, z_n, mapdata] = loewner_slider(direction, tooth, z_n, w_n, mapdata)

persistent unzip zipup compute_a
if isempty(unzip)
  from shapelab.loewner import compute_driving_step as unzip
  from shapelab.loewner import compute_a
end

switch direction
case 'down'
  if tooth==1 % Then initialize data array
    N = length(z_n) - 3;
    mapdata.a_array = 0;  % Allocation for the terminal Moebius map
    mapdata.loewner_data.lambda = zeros([(N-2) mapdata.M+1]);
    mapdata.loewner_data.dlambda = zeros([(N-2) mapdata.M]);
    mapdata.loewner_data.s = zeros([(N-2) mapdata.M]);
    lambda = 0;
    s = 0;
  else
    lambda = mapdata.loewner_data.lambda(tooth-1, end);
    s = mapdata.loewner_data.s(tooth-1, end);
  end

  a_id = tooth + 2;
  a = z_n(a_id);

  %gn = [w_n(1:a_id-2); z_n(a_id+1:end)];
  gn = w_n;
  %g = [z_n(1:a_id-2); z_n(a_id+1:end)];
  g = z_n;

  a = compute_a(g, lambda);
  an = compute_a(gn, lambda);

  % Now unzip the slit
  [lambdan, dlambdan, ds, a, an, g, gn] = unzip(a, an, g, gn, lambda, a_id-1, mapdata.M);

  mapdata.loewner_data.lambda(tooth, :) = lambdan;
  mapdata.loewner_data.dlambda(tooth, :) = dlambdan;
  mapdata.loewner_data.s(tooth, :) = s + cumsum(ds);

 % interior pts , 'left' points , 'right' points
  %z_n(1:a_id-2) = g(1:a_id-2);
  %w_n(1:a_id-2) = gn(1:a_id-2);
  z_n = g;
  w_n = gn;

  %z_n(a_id+1:end) = g(a_id-1:end);
  %w_n(a_id+1:end) = gn(a_id-1:end);

  %z_n(a_id) = 0;
  %w_n(a_id) = 0;

case 'up'
  error('Not yet implemented');
end
