function[] = visualize(self, z, z_in, z_out, visdata)
% visualize -- visualization for zipper drivers

%persistent moebius moebius_plot
%
%if isempty(moebius)
%  from shapelab.common import moebius
%  moebius_plot = [-1, i;...
%                   1, i];  % Just a map from half plane to unit circle.
%end

if visdata.visualize
  %z = moebius([z],moebius_plot);
  z = self.moebius_maps.H_to_D(z);
  %z_in = moebius(z_in,moebius_plot);
  z_in = self.moebius_maps.H_to_D(z_in);
  %z_out = moebius(z_out,moebius_plot);
  z_out = self.moebius_maps.H_to_D(z_out);

  set(visdata.shape_plot, 'xdata', real(z), 'ydata', imag(z));
  set(visdata.zin_plot, 'xdata', real(z_in), 'ydata', imag(z_in));
  set(visdata.zout_plot, 'xdata', real(z_out), 'ydata', imag(z_out));

  drawnow;
end
