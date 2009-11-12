function[z]= evaluate_inverse_map(w,mapdata,varargin)
% [z]= evaluate_inverse_map(w,mapdata,{point_id=zeros(size(w))})
%
%     Evaluates the inverse of the zipper-type conformal map defined by mapdata.
%     The optional input 'point_id' should be specified as follows:
%
%     z(point_id==0) -----> (default) points somewhere inside or outside the
%     unit circle that you want to map inside or outside the shape.
%
%     z(point_id==1) -----> points on the interior boundary of the unit disc
%     that you want to map onto the shape. 
%
%     z(point_id==2) -----> points on the exterior boundary of the unit disc
%     that you want to map onto the shape.
%
%     [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for
%          conformal mapping", 2006.

persistent input_schema zip moebius moebius_inv csqrt pcpow ncpow
if isempty(input_schema)
  from labtools import input_schema
  from shapelab.common import moebius
  from shapelab.common import moebius_inverse as moebius_inv
  from shapelab.common import positive_angle_square_root as csqrt
  from shapelab.common import positive_angle_exponential as pcpow
  from shapelab.common import negative_angle_exponential as ncpow

  imp shapelab.zipper as zip
end

opt = input_schema({'point_id'}, {zeros(size(w))}, [], varargin{:});
switch lower(mapdata.type)
case 'geodesic'
  ifa = zip.geodesic.inverse_base_conformal_map;
  zipper = false;
case 'slit'
  ifa = zip.slit.inverse_base_conformal_map;
  zipper = false;
case 'zipper'
  warning('For exterior points, I can''t guarantee that this won''t return garbage');
  ifa = zip.zipper.inverse_base_conformal_map;
  zipper = true;
case 'zipper_weld'
  ifa = zip.zipper.inverse_base_conformal_map;
  ifa_geo = zip.geodesic.inverse_base_conformal_map;
  zipper = true;
end


% We must deal with three cases:
w0 = w(opt.point_id==0);
w1 = w(opt.point_id==1);
w2 = w(opt.point_id==2);

z = zeros(size(w));

%% Hokay, here we go

% If we're on the interior/exterior proper, no hocus-pocus is necessary
if any(opt.point_id==0)
  % Meh, just in case
  if any(abs(w0)==1)
    error('Point(s) specified as strictly not on unit circle looks like it is');
  end
  ext_points = abs(w0)>1;
  int_points = abs(w0)<1;

  % map unit circle to R:
  m = [-1, i;...
        1, i];
  w0 = moebius_inv(w0,m);

  % apply inv(m_in) on interior
  w0(int_points) = moebius_inv(w0(int_points), mapdata.m_in);

  % apply inv(m_out) on exterior
  w0(ext_points) = moebius_inv(w0(ext_points), mapdata.m_out);

  % Invert the terminal map
  m = [1, 0;...
       -mapdata.a_array(end), 1];
  if strcmpi(mapdata.type, 'zipper')
    alpha = pi - angle(moebius(mapdata.c_array(end),m));
    %w0 = exp(i*(pi-alpha))*(w0*sign(mapdata.winding_number)).^(alpha/pi);
    w0 = w0*sign(mapdata.winding_number);
    % regrettably, the non-conformal terminal map requires different treatment
    % of interior and exterior points
    w0(int_points) = exp(i*(pi-alpha))*pcpow(w0(int_points), 'alpha', alpha/pi);
    w0(ext_points) = exp(i*(pi-alpha))*ncpow(w0(ext_points), 'alpha', alpha/pi);
    w0 = moebius_inv(w0,m);
  elseif strcmpi(mapdata.type, 'zipper_weld')
    % Do the geodesic/slit terminal map inverse + 1 inverse geodesic base map
    w0 = csqrt(-sign(mapdata.winding_number)*w0);
    w0 = moebius_inv(w0,m);

    % 1 geodesic inverse base map
    ifa_opt.point_id = zeros(size(w0));
    ifa_opt.cut_magnitude = mapdata.zip_magnitude;
    w0 = ifa_geo(w0, mapdata.c_array(end), ifa_opt);
  else
    % This should map H (interior) ----> left-half plane (right-half if winding
    % number is negative). Technically, this is a special case of the zipper
    % type. Meh.
    w0 = csqrt(-sign(mapdata.winding_number)*w0);
    w0 = moebius_inv(w0,m);
  end

  % Now loop through a_array
  N = length(mapdata.a_array)+1;

  ifa_opt.point_id = zeros(size(w0));
  ifa_opt.cut_magnitude = mapdata.zip_magnitude;
  for q = mapdata.N_teeth:-1:1
    if zipper
      w0 = ifa(w0, mapdata.c_array(q), mapdata.a_array(q), ifa_opt);
    else
      w0 = ifa(w0, mapdata.a_array(q), ifa_opt);
    end
  end

  % Now it's easy-peasy:
  if zipper
    w0 = moebius_inv(w0.^2, mapdata.m_initial);
  else
    w0 = moebius_inv((w0/i).^2, mapdata.m_initial);
  end
  z(opt.point_id==0) = w0;
end

% If we're on the unit circle interior
if any(opt.point_id==1)
  if not(all(abs(abs(w1)-1)<1e-10))
    error('It doesn''t look like you gave me points on the unit circle....');
  end
  % map unit circle to R:
  m = [-1, i;...
        1, i];
  w1 = real(moebius_inv(w1,m));

  % apply inv(m_in) on interior
  w1 = moebius_inv(w1, mapdata.m_in);

  if strcmpi(mapdata.type, 'zipper')
    m = [1 0; ...
         -mapdata.a_array(end) 1];
    alpha = pi - angle(moebius(mapdata.c_array(end),m));
    w1 = w1*sign(mapdata.winding_number);
    temp = w1;
    w1(temp<0) = -abs(w1(temp<0)).^(alpha/pi);
    w1(temp>=0) = exp(i*(pi-alpha))*(w1(temp>=0)).^(alpha/pi);
    w1 = moebius_inv(w1,m);
  elseif strcmpi(mapdata.type, 'zipper_weld')
    % Do the geodesic/slit terminal map inverse + 1 inverse geodesic base map
    w1 = w1*-sign(mapdata.winding_number);
    temp = w1;
    w1(temp<0) = i*sqrt(abs(temp(temp<0)));
    w1(temp>=0) = -sqrt(temp(temp>=0));
    % Invert the terminal map (moebius)
    m = [1, 0;...
         -mapdata.a_array(end), 1];
    w1 = moebius_inv(w1,m);

    % 1 geodesic inverse base map
    ifa_opt.point_id = ones(size(w1));
    ifa_opt.point_id(abs(imag(w1))>1e-10) = 0;
    ifa_opt.cut_magnitude = mapdata.zip_magnitude;
    w1 = ifa_geo(w1, mapdata.c_array(end), ifa_opt);
  else
    w1 = w1*-sign(mapdata.winding_number);
    temp = w1;
    w1(temp<0) = i*sqrt(abs(temp(temp<0)));
    w1(temp>=0) = -sqrt(temp(temp>=0));
    % Invert the terminal map (moebius)
    m = [1, 0;...
         -mapdata.a_array(end), 1];
    w1 = moebius_inv(w1,m);
  end

  ifa_opt.point_id = ones(size(w1));
  ifa_opt.point_id(abs(imag(w1))>1e-10) = 0;
  ifa_opt.cut_magnitude = mapdata.zip_magnitude;

  % Here's the fun part: right now everything's on the real line. The second it
  % gets zipped up from the real line, we have to tell ifa to treat it as an
  % interior point instead of a boundary point.

  N = length(mapdata.a_array)+1;
  for q = (N-2):-1:1
    if zipper
      w1 = ifa(w1, mapdata.c_array(q), mapdata.a_array(q), ifa_opt);
    else
      w1 = ifa(w1, mapdata.a_array(q), ifa_opt);
    end
    winterior = abs(imag(w1))>0;
    ifa_opt.point_id(winterior) = 0;
  end

  if zipper
    w1 = moebius_inv(w1.^2, mapdata.m_initial);
  else
    w1 = moebius_inv((w1/i).^2, mapdata.m_initial);
  end
  z(opt.point_id==1) = w1;
end

% If I cared more, I could do this without copy+paste from above -- only two
% lines are different
% If we're on the unit circle exterior
if any(opt.point_id==2)
  if not(all(abs(abs(w2)-1)<1e-10))
    error('It doesn''t look like you gave me points on the unit circle....');
  end
  % map unit circle to R:
  m = [-1, i;...
        1, i];
  w2 = real(moebius_inv(w2,m));

  % apply inv(m_out) on exterior
  w2 = moebius_inv(w2, mapdata.m_out);

  if strcmpi(mapdata.type, 'zipper')
    m = [1 0; ...
         -mapdata.a_array(end) 1];
    alpha = pi - angle(moebius(mapdata.c_array(end),m));
    w2 = w2*sign(mapdata.winding_number);
    temp = w2;
    w2(temp<0) = exp(i*(pi-alpha))*(w2(temp<0)).^(alpha/pi);
    w2(temp>=0) = w2(temp>=0).^(alpha/pi);
    %w2(temp>=0) = sqrt(w2(temp>=0));
    w2 = moebius_inv(w2,m);
  elseif strcmpi(mapdata.type, 'zipper_weld')
    % Do the geodesic/slit terminal map inverse + 1 inverse geodesic base map
    w2 = w2*-sign(mapdata.winding_number);
    temp = w2;
    w2(temp<0) = i*sqrt(abs(w2(temp<0)));
    w2(temp>=0) = sqrt(w2(temp>=0));
    % Invert the terminal map (moebius)
    m = [1, 0;...
         -mapdata.a_array(end), 1];
    w2 = moebius_inv(w2,m);

    % 1 geodesic inverse base map
    ifa_opt.point_id = ones(size(w2));
    ifa_opt.point_id(abs(imag(w2))>1e-10) = 0;
    ifa_opt.cut_magnitude = mapdata.zip_magnitude;
    w2 = ifa_geo(w2, mapdata.c_array(end), ifa_opt);
  else
    w2 = w2*-sign(mapdata.winding_number);
    temp = w2;
    w2(temp<0) = i*sqrt(abs(w2(temp<0)));
    w2(temp>=0) = sqrt(w2(temp>=0));
    % Invert the terminal map (moebius)
    m = [1, 0;...
         -mapdata.a_array(end), 1];
    w2 = moebius_inv(w2,m);
  end



  ifa_opt.point_id = ones(size(w2));
  ifa_opt.point_id(abs(imag(w2))>1e-10) = 0;
  ifa_opt.cut_magnitude = mapdata.zip_magnitude;

  % Here's the fun part: right now everything's on the real line. The second it
  % gets zipped up from the real line, we have to tell ifa to treat it as an
  % interior point instead of a boundary point.

  N = length(mapdata.a_array)+1;
  for q = (N-2):-1:1
    if zipper
      w2 = ifa(w2, mapdata.c_array(q), mapdata.a_array(q), ifa_opt);
    else
      w2 = ifa(w2, mapdata.a_array(q), ifa_opt);
    end
    winterior = abs(imag(w2))>0;
    ifa_opt.point_id(winterior) = 0;
  end

  if zipper
    w2 = moebius_inv(w2.^2, mapdata.m_initial);
  else
    w2 = moebius_inv((w2/i).^2, mapdata.m_initial);
  end
  z(opt.point_id==2) = w2;
end
