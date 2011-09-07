function[visdata]=initialize_visualization(self, vistf, z)
% initialize_visualization -- Initializes visualization parameters
%
% visdata = initialize_visualization(self, vistf, z)

if vistf
  figure(); 
  theta = linspace(0,2*pi,200);
  plot(exp(i*theta), 'b--');
  axis square;
  axis off;
  hold on;
  H = self.moebius_maps.H_to_D;
  %shape_plot = plot(moebius(z, moebius_maps.H_to_D), 'k.-');
  shape_plot = plot(H(z), 'k.-');
  %temp = moebius([Inf; 0],moebius_maps.H_to_D);
  temp = H([Inf; 0]);
  zin_plot = plot(real(temp), imag(temp), 'r.');
  %temp = moebius([Inf; 0], moebius_maps.H_to_D);
  temp = H([Inf; 0]);
  zout_plot = plot(real(temp), imag(temp), 'b.');
else
  shape_plot = 0;
  zin_plot = 0;
  zout_plot = 0;
end
                                                               
visdata.shape_plot = shape_plot;
visdata.zin_plot = zin_plot;
visdata.zout_plot = zout_plot;
