function[H] = ba_fingerprint_driver(stuff, z)
% ba_fingerprint -- Evaluates the Beurling-Ahlfors extension map
%
% H = ba_fingerprint_driver(stuff, z)
%
%     Given fingerprint information in the struct stuff, this evaluates the
%     BA extension ***ON H*** at the locations z (also on H). See
%     ba_fingerprint and weld_primitive_setup for details.

persistent wextend
if isempty(wextend)
  from shapelab.extensions import weld_primitive_extend_driver as wextend
end

% H = 1/(2*y)*(a*Iplus + conj(a)*Iminus)
a = 1+2*i;
abar = conj(a);

x = real(z);
y = imag(z);

conj_flags = y<0;
y(conj_flags) = -y(conj_flags);

Ix = wextend(x, stuff);
Iplus = wextend(x+y, stuff) - Ix;
Iminus = Ix - wextend(x-y, stuff);

H = 1./(2*y).*(a*Iplus + abar*Iminus);
H(conj_flags) = conj(H(conj_flags));
