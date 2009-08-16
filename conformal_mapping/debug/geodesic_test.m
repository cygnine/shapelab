global handles;
gd = handles.shapelab.conformal_mapping.geodesic;

N = 100;
theta = linspace(-pi,pi,N+1); 
theta = theta(1:N);
z = (1+0.75*cos(theta*4)).*exp(i*theta);

thetaf = linspace(-pi,pi,10*N);
zf = (1+0.75*cos(thetaf*4)).*exp(i*thetaf);
M = 50;
wint = zeros([N*M,1]);
for q = 1:M
  wint((q-1)*length(z)+1:q*length(z)) = z.*((q)/(M+1));
end

wout = zeros([N*M,1]);
for q = 1:M
  wout((q-1)*length(z)+1:q*length(z)) = z.*((q+M+1)/M); 
end

mapin = gd.compute_map_coordinates(z,'z_in', 0, 'z_out', 0);

%mapout = gd.compute_map_coordinates(fliplr(z), 'z_in', Inf, 'z_out', Inf);

fwint = gd.evaluate_map(wint,mapin);
%fwout = gd.evaluate_map(wout,mapout);

tin = gd.evaluate_map(z,mapin);
wshape = gd.evaluate_inverse_map(exp(i*thetaf), mapin);
tinf = gd.evaluate_map(zf,mapin);
%tout = gd.evaluate_map(z,mapout);

%%% Debugging wtf the csqrt mapping doesn't do what it's supposed to
mapdata = mapin;

fa = handles.shapelab.conformal_mapping.geodesic.base_conformal_map;
moebius = handles.shapelab.common.moebius;
dab = handles.shapelab.common.disc_a_to_b;

[z_initial, a_array, a_cut_bias, normalization_persistence, zeta_n, normalization] = ...
  deal(mapdata.z_initial, mapdata.a_array, mapdata.a_cut_bias, ...
    mapdata.normalization_persistence, mapdata.zeta_n, mapdata.normalization);

% Initial map
w = i*sqrt(moebius(wint,[1 -z_initial(2); 1 -z_initial(1)]));
w(not(isfinite(w))) = Inf;
figure; p1 = plot(w, 'b.');

for q = 1:length(a_array)
  w = fa(w,a_array(q));
  if q<length(a_array) && normalization_persistence(q)~=0
    w = w*normalization_persistence(q);
  end
  set(p1, 'xdata', real(w), 'ydata', imag(w));
  drawnow; pause;
end

% Terminal map:
w = -(moebius(w,[1 0; -1/zeta_n 1])).^2;

% Map to unit circle:
w = moebius(w, [-1 i; 1 i]);

% Normalize
w = dab(w,normalization(1), normalization(2));
