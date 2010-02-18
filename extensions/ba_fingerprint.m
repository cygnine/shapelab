function[H] = ba_fingerprint(theta_int, theta_ext, z)
% ba_fingerprint -- Evaluates the Beurling-Ahlfors extension map
%
% H = ba_fingerprint(theta_int, theta_ext, z)
%
%     Assuming that the fingerprint is periodically extended over R, this
%     function takes input values (theta_int, theta_ext) of the fingerprint and
%     evaluates the Beurling-Ahlfors extension of the fingerprint h at the
%     locations z. Any inputs z that do not lie in the upper half-plane are
%     evaluated as: H <----- conj(H(conj(z)))
%
%     The underlying approximation is a piecewise polynomial -- see
%     weld_primitive_extend. The input fingerprint must satisfy theta_int(1) =
%     theta_ext(1) = 0 and max(theta)<(2*pi).

persistent wsetup jacobi_glq ba
if isempty(wsetup)
  from shapelab.extensions import weld_primitive_setup as wsetup
  from speclab.orthopoly1d.jacobi.quad import gauss_lobatto_quadrature as jacobi_glq
  from shapelab.extensions import ba_fingerprint_driver as ba
end

stuff = wsetup(theta_int, theta_ext, 'local_nodes', jacobi_glq(9));
H = ba(stuff, z);
