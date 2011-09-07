function[w, z, self] = unzip(self, tooth, w, z)
% unzip -- Unzips one tooth using the geodesic algorithm
%
% [w, z, self] = unzip(self, tooth, w, z)
%
%     w -- points on the interior
%     z -- points on the exterior

persistent zipdown
if isempty(zipdown)
  from shapelab.loewner.solutions import normal_linear_slit_unzip as zipdown
end

%if tooth==1 % Then initialize data array
%  N = length(z_n) - 3;
%  mapdata.a_array = zeros([N-1 1]);
%end

a_id = tooth + 2;
a = z(a_id);
self.a_array(tooth) = a;

b = abs(a)^2/real(a);  % The point on real axis that the next map sends to infinity
c = abs(a)^2/imag(a);  % The image of a after the map

% Use a moebius map m to send the circular arc 0 -- a to the perpendicular line 
% 0 -- mapdata.tooth_length/c
map = MoebiusMap([self.tooth_length/c 0; -1/b 1]);
self.moebius_maps.tooth_maps{tooth} = map;
zinf = z(end-2);
z = map(z);
w = map(w);

% Update deriative at infinity
self.derivative_at_inf = self.derivative_at_inf*map.derivative(zinf);
%H = map.H;
%self.derivative_at_inf = self.derivative_at_inf*...
%    det(H)/(H(2,1)*zinf + H(2,2))^2;

% Now unzip the perpendicular line segment
[temp, temp2, temp3] = zipdown(i*self.tooth_length, ...
                               [z(1:a_id-2); w(1:a_id-2); z(a_id+1:end)], ...
                               [z(a_id-1:a_id); w(a_id-1:a_id)]);

zinf = z(end-2);

% interior pts , 'left' points , 'right' points
z(1:a_id-2) = temp(1:a_id-2);
w(1:a_id-2) = temp(a_id-1:2*a_id-4);

z(a_id+1:end) = temp(2*a_id-3:end);
w(a_id+1:end) = temp(2*a_id-3:end);

z(a_id-1:a_id) = temp3(1:2);
w(a_id-1:a_id) = temp2(3:4);

z(a_id) = 0;
w(a_id) = 0;

% update derivative
%self.derivative_at_inf = self.derivative_at_inf*...
%     (zinf-real(a))/(z(end-2) - real(a));
self.derivative_at_inf = self.derivative_at_inf*...
     (zinf)/(z(end-2));
