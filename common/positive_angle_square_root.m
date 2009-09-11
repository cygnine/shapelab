function[w] = positive_angle_square_root(z,varargin)
% [w] = positive_angle_square_root(z,{cut_bias=[]})
% 
%     The input z (vector-supported) has form r*exp(i*phi), where phi \in
%     [0,2*pi). This function returns sqrt(r)*exp(i*phi/2). Note that this
%     produces a branch discontinuity at phi = 2*pi. I.e., the branch cut of the
%     sqrt function is placed on the positive real axis.
%
%     The optional input cut_bias specifies which side of the branch the
%     positive x-axis gets mapped to, and if this option is set then some small
%     tolerance leeway is given in favor of the bias. cut_bias=true means the
%     positive x-axis gets mapped to itself. cut_bais=false means the positive
%     x-axis gets mapped to the negative x-axis.
%
%     If cut_bias is an array, it is interpreted as a vector of boolean
%     (true=positive, false=negative) indicators determining which inputs take
%     the (+) x-axis side of the branch and which take the (-) x-axis side of
%     the branch. For inputs not on the x-axis, this indicator is irrelevant.
%
%     In light of positive_angle_exponential, this function will probably be
%     deprecated soon.

global handles;
opt = handles.common.input_schema({'cut_bias'}, {[]}, [], varargin{:});

if length(opt.cut_bias)==1
  if opt.cut_bias
    tol = 1e-8;
  else
    tol = -1e-8;
  end
else
  zflags = not(opt.cut_bias);
  tol = 0;
end

ang = mod(angle(z),2*pi);
flags = ang<-tol;
ang(flags) = ang(flags)+2*pi;
w = sqrt(abs(z)).*exp(i*ang/2);
