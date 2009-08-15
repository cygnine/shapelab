function[w] = positive_angle_square_root(z)
% [w] = positive_angle_square_root(z)
% 
%     The input z (vector-supported) has form r*exp(i*phi), where phi \in
%     [0,2*pi). This function returns sqrt(r)*exp(i*phi/2). Note that this
%     produces a branch discontinuity at phi = 2*pi.

ang = angle(z);
flags = ang<0;
ang(flags) = ang(flags)+2*pi;
w = sqrt(abs(z)).*exp(i*ang/2);
