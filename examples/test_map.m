global handles;
gd = handles.ConformalMapping.GeodesicAlgorithm;
explot = handles.common.explot;

N = 100;
theta = linspace(0,2*pi,N+1);
theta = theta(1:N).';
z = exp(i*theta);

%[z_initial,a_array,zeta_n] = gd.compute_map_coordinates(z);

[z_initial2,a_array2,zeta_n2] = gd.visualize_map_computation(z);

w = gd.evaluate_map(z,z_initial,a_array,zeta_n);
