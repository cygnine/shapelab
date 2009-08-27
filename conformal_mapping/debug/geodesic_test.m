global handles;
gd = handles.shapelab.conformal_mapping.geodesic;

N = 100;
theta = linspace(-pi,pi,N+1); 
theta = theta(1:N);
z = (1+0.45*cos(theta*4)).*exp(i*theta);

thetaf = linspace(-pi,pi,10*N);
zf = (1+0.45*cos(thetaf*4)).*exp(i*thetaf);
M = 50;
zint = zeros([N*M,1]);
for q = 1:M
  zint((q-1)*length(z)+1:q*length(z)) = z.*((q)/(M+1));
end

zout = zeros([N*M,1]);
for q = 1:M
  zout((q-1)*length(z)+1:q*length(z)) = z.*((q+M+1)/M); 
end

mapdata = gd.compute_map_coordinates(z,'z_in', 0, 'w_in', 0, ...
      'zip_magnitude', 0.5);

wint = zeros([N*M 1]);
for q = 1:M
  wint((q-1)*length(z)+1:q*length(z)) = exp(i*theta).*(q-1)/M;
end

wout = zeros([N*M 1]);
for q = 1:M
  wout((q-1)*length(z)+1:q*length(z)) = exp(i*theta).*(q+M)/M;
end

wint_image = gd.evaluate_inverse_map(wint,mapdata);
wout_image = gd.evaluate_inverse_map(wout,mapdata);

% testing:
%test1 = gd.evaluate_map(z,mapin);
%test2 = gd.evaluate_map(z,mapout);
%
%test3 = gd.evaluate_inverse_map(test1,mapin);
%test4 = gd.evaluate_inverse_map(test2,mapout,'boundary', true);
%
%fwint = gd.evaluate_map(wint,mapin);
%fwout = gd.evaluate_map(wout,mapout);

%tin = gd.evaluate_map(z,mapin);
%wshape = gd.evaluate_inverse_map(exp(i*thetaf), mapin);
%tinf = gd.evaluate_map(zf,mapin);
%tout = gd.evaluate_map(z,mapout);

% Debugging stuff
moebius = handles.shapelab.common.moebius;
fa = handles.shapelab.conformal_mapping.geodesic.base_conformal_map;
opt.z_in = 0;
opt.z_out = 0;
opt.winding_number = 1;

z_n = z(:);
N = length(z_n);
assert(N>2,'Error: you must give at least three points defining a shape');
% Let's throw infinity in there for kicks
z_n(end+1) = Inf;
wint = zint(:); % images of interior point
z_n(end+1) = opt.z_in;

%% These are needed as output for the initial map
z_initial = [z_n(1); z_n(2)];

%% Temporary array, order is now z(2),z(3),...,z(n),Inf,z_norm,z(1)
zeta_n = circshift(z_n,-1);
% This is weird to me: taking the sqrt branch cut at pi instead of 2*pi. I don't
% understand...just following Marshall. If I use a sqrt branch cut at 2*pi (as
% is done in the rest of the algorithm), I shouldn't need the i factor.
m = [1 -z_initial(2);...
     1 -z_initial(1)];
zeta_n = i*sqrt(moebius(zeta_n(2:end), m)); % don't need original z(2)
zeta_n(end) = Inf; % sadly, sqrt(complex Inf) = NaN
wint = i*sqrt(moebius(wint, m));

unzipped_in = [Inf;0]; % original z_1, z_2
unzipped_out = [Inf;0];
%figure; 
%p1 = plot(zeta_n, 'b.'); hold on;
%p2 = plot(wint, 'g.');
%p3 = plot(unzipped_in, 'r.');
%p4 = plot(unzipped_out, 'k.');

% Ok, now zip down each point, keeping track of where everything goes
for q = 1:(N-2);
  % The next map is defined by the parameter:
  a_array(q) = zeta_n(1);

  %% Apply the map to all the points
  % interior points + infinity + normalization point + original z(1)
  temp = zeros(size(zeta_n(2:end)));
  temp(end) = 2;
  zeta_n = fa(zeta_n(2:end),a_array(q),'point_id',temp,'cut_magnitude', 1/2);
  %zeta_n = fa(zeta_n(2:end),a_array(q),'point_id',temp);
  wint = fa(wint,a_array(q),'point_id',zeros(size(wint)),'cut_magnitude', 1/2);
  %wint = fa(wint,a_array(q),'point_id',zeros(size(wint)));
  temp = 2*ones(size(unzipped_in), 'int8'); temp(end) = 1;

  unzipped_in_test1 = unzipped_in;

  [unzipped_in,blah] = fa(unzipped_in,a_array(q),'point_id',temp,'cut_magnitude',1/2); 
  %[unzipped_in,blah] = fa(unzipped_in,a_array(q),'point_id',temp);
  [blah,unzipped_out] = fa(unzipped_out,a_array(q),'point_id',temp,'cut_magnitude',1/2);
  %[blah,unzipped_out] = fa(unzipped_out,a_array(q),'point_id',temp);

  % Now add 0 to the list of zipped_in/out points:
  unzipped_in(end+1) = 0;
  unzipped_out(end+1) = 0;

  %set(p1, 'xdata', real(zeta_n), 'ydata', imag(zeta_n));
  %set(p2, 'xdata', real(wint), 'ydata', imag(wint));
  %set(p3, 'xdata', real(unzipped_in), 'ydata', imag(unzipped_in));
  %set(p4, 'xdata', real(unzipped_out), 'ydata', imag(unzipped_out));
  %drawnow;
end

unzipped_in_test2 = unzipped_in;

% Finally, the terminal map. By now, only 3 points are left: infinity +
% normalization point + z(1)

a_array(end+1) = 1/zeta_n(3);
m = [ 1              0; ...
     -a_array(end) 1];
%% Care about where infinity + normalization point go:
% Infinity is on the outside
zeta_n(1) = -sign(opt.winding_number)*moebius(zeta_n(1), m).^2;
% normalization point is on the inside
zeta_n(2) = -sign(opt.winding_number)*moebius(zeta_n(2), m).^2;
% map zipped points as well:
unzipped_in = -sign(opt.winding_number)*moebius(unzipped_in,m).^2;
unzipped_out = -sign(opt.winding_number)*moebius(unzipped_out,m).^2;

% Now zeta_n(1) is probably not infinity...let's force it to be:
m1 = [-1, i; ...
       1, i];  % to unit circle
m2 = [0, 1;...
      1 0]; % inversion so unit disc == outside
a = moebius(zeta_n(1), m2*m1); % image of infinity to disc
m3 = [abs(a)/a -abs(a); ...
      -conj(a) 1];  % map a (original infinity) to 0
m4 = [0, -1;...
      -1 0]; % invert m2
m5 = [i, -i; ...
      -1, -1]; % invert m1
% total map:
m_out = m5*m4*m3*m2*m1;
m_out = real(m_out/m_out(1,1));

m6 = m1;
a = moebius(zeta_n(2),m1);
m7 = [abs(a)/a -abs(a); ...
      -conj(a) 1];  % map normalization point to 0
m8 = m5;
m_in = m8*m7*m6;
m_in = real(m_in/m_in(1,1));

% Total map
unzipped_in = moebius(unzipped_in,m_in);
unzipped_out = moebius(unzipped_out,m_out);

% Finally, map to unit circle
m = [-1, i;...
      1, i];
unzipped_in = moebius(unzipped_in, m);
unzipped_out = moebius(unzipped_out, m);

tin = unwrap(angle(unzipped_in));
tout = unwrap(angle(unzipped_out));

wincircle_image = gd.evaluate_inverse_map(unzipped_in, mapdata,...
    'point_id', ones(size(unzipped_in)));

woutcircle_image = gd.evaluate_inverse_map(unzipped_out, mapdata,...
    'point_id', 2*ones(size(unzipped_out)));
