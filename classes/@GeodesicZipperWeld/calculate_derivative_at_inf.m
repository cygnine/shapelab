function[val] = calculate_derivative_at_inf(self, zinf)
% calculate_derivative_at_inf -- Calculates derivative of the map at infinity
%
% val = calculate_derivative_at_inf(self)
%
%     Using the map defined by self, this computes the derivative at
%     infinity.

persistent zipdown
if isempty(zipdown)
  from shapelab.loewner.solutions import normal_linear_slit_unzip as zipdown
end

if nargin<2
  zinf = Inf;
end
flags = isinf(zinf);

% Initial map
H = self.moebius_maps.initial_map.H;
val(flags) = det(H)/H(2,1)^2;
val(not(flags)) = det(H)./(H(2,1)*zinf(not(flags)) + H(2,2)).^2;

zinf(flags) = H(1,1)/H(2,1);
zinf(not(flags)) = self.moebius_maps.initial_map(zinf(not(flags)));

val = val*i/2./sqrt(zinf);
zinf = i*sqrt(zinf);

% Tooth maps
for q = 1:self.N_teeth
  H = self.moebius_maps.tooth_maps{q}.H;

  val = val.*det(H)./(H(2,1)*zinf + H(2,2)).^2;
  zinf = self.moebius_maps.tooth_maps{q}(zinf);

  zback = zinf;
  zinf = zipdown(i*self.tooth_length, zinf);
  val = val.*zback./zinf;
end

% Terminal map
H = self.moebius_maps.terminal_map.H;
val = val.*det(H)./(H(2,1)*zinf + H(2,2)).^2;
zinf = self.moebius_maps.terminal_map(zinf);

val = val.*-2*sign(self.winding_number).*zinf;
zinf = -sign(self.winding_number).*zinf.^2;

% Moebius alignment
H = self.moebius_maps.exterior_terminal.H;
val = val.*det(H)./(H(2,1)*zinf + H(2,2)).^2;
zinf = self.moebius_maps.exterior_terminal(zinf);

% Now zinf better be at -i
asfd

% Final H_to_D factor
H = self.moebius_maps.H_to_D.H;
val(flags) = val(flags)*(-2*i);
val(not(flags)) = val(not(flags)).*det(H)./(H(2,1)*zinf(not(flags)) + H(2,2)).^2;

val = det(self.moebius_maps.exterior_rotation.H)*val;
