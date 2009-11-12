function[z_initial,a_array,zeta_n] = visualize_map_computation(z_n)
% [Z_INITIAL,A_ARRAY,ZETA_N] = VISUALIZE_MAP_COMPUTATION(Z_N)
%
%     Computes an array of complex values A_ARRAY, with length N-2 where N is
%     the length of the input vector Z_N. Also returns Z_INITIAL=[Z_N(1),
%     Z_N(2)], and ZETA_N, the location of the final node. These data
%     corresponds to the parameters defining the composition maps of the
%     `geodesic' algorithm for conformal mapping. See [1]. This code is
%     identical to COMPUTE_MAP_COORDINATES except that it visualizes the map
%     formation; the style for this is shamelessly stolen from M. Feiszli.
%
%     The initial mapping sending Z_N(1) off to infinity, and the terminal
%     mapping sending the real line to the unit circle are assumed to be
%     understood without reference. See implementation of EVALUATE_MAP.
%
%  [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%  mapping", 2006.

error('This function seems to be deprecated');
global packages;
evaluate_fa = packages.ConformalMapping.GeodesicAlgorithm.evaluate_fa;
figure();
hold on;

z_n = z_n(:);
N = length(z_n);
assert(N>2,'Error: you must give at least three points defining a shape');

a_array = zeros([N-3,1]);

myplot = plot(z_n, 'r.');
axis([-1,1,-1,1]);
axis off;

% These are needed as output for the initial map
z_initial = [z_n(1); z_n(2)];

% Temporary array, first send z_n(1) to infinity
zeta_n = circshift(z_n,-1);
zeta_n = i*sqrt((zeta_n - z_n(2))./(zeta_n - z_n(1)));

for q = 1:(N-2)
  % Determine next parameter a
  a_array(q) = zeta_n(q+1);

  % Operate on remaining points with f_a
  zeta_n = evaluate_fa(zeta_n,a_array(q));
  
  % Turns out Matlab doesn't like Infs
  if q==1
    zeta_n(end) = -abs(a_array(1))^2/real(a_array(1));
  end
  plot_zeta = (i-zeta_n)./(i+zeta_n);
  set(myplot,'xdata',real(plot_zeta));
  set(myplot,'ydata',imag(plot_zeta));
  axis([-1,1,-1,1]);
  drawnow();
end

% Terminal map:
plot_zeta = -(zeta_n./(1-zeta_n./zeta_n(end))).^2;
plot_zeta = (i-plot_zeta)./(i+plot_zeta);
plot_zeta(end) = -1;
set(myplot,'xdata',real(plot_zeta));
set(myplot,'ydata',imag(plot_zeta));
drawnow();

zeta_n = zeta_n(end);
