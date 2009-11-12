function[w] = negative_angle_exponential(z,varargin)
% [w] = negative_angle_exponential(z,{alpha=1/2, cut_bias=[]})
% 
%     The input z (vector-supported) has form r*exp(i*-phi), where phi \in
%     [0,2*pi). This function returns r^(alpha)*exp(i*-phi*alpha). Note that this
%     produces a branch discontinuity at phi = 2*pi. I.e., the branch cut of the
%     power function is placed on the positive real axis. However, in contrast
%     to positive_angle_exponential, the branch cut happens from angle(z) =
%     -2*pi+eps to angle(z) = -2*pi - eps.
%
%     The optional input cut_bias specifies which side of the branch the
%     positive x-axis gets mapped to, and if this option is set then some small
%     tolerance leeway is given in favor of the bias. cut_bias=true means the
%     positive x-axis gets mapped to itself. cut_bais=false means the positive
%     x-axis gets mapped to the line segment z=exp(i*2*pi*alpha).
%
%     If cut_bias is an array, it is interpreted as a vector of boolean
%     (true=positive, false=negative) indicators determining which inputs take
%     the (+) x-axis side of the branch and which take the (-) x-axis side of
%     the branch. For inputs not on the x-axis, this indicator is irrelevant.

persistent input_schema
if isempty(input_schema)
  from labtools import input_schema
end

opt = input_schema({'cut_bias', 'alpha'}, {true,1/2}, [], varargin{:});

if length(opt.cut_bias)==1
  if opt.cut_bias
    tol = 1e-14;
  else
    tol = -1e-14;
  end
else
  zflags = not(opt.cut_bias);
  tol = 0;
end

ang = mod(angle(z),-2*pi);
flags = ang<(-2*pi+tol);
ang(flags) = ang(flags)+2*pi;
w = (abs(z)).^opt.alpha.*exp(i*ang*opt.alpha);
