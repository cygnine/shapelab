function g = invert_a(a, xi, branch);
% invert_a -- inverts the g(a) relation
%
% g = invert_a(a, xi, branch)
%
%     The relation a = g^2 - 2*g*xi is inverted to solve for g, with appropriate
%     choices of branches for g. The branch in the upper-half plane is always
%     chosen.
%
%     When a is complex, the branch is always chosen automatically; when a is
%     real, the boolean array branch determines which solution is chosen. true
%     corresponds to the greater root, and false the less root.

reals = (imag(a)==0);

discriminant = sqrt(2*a + xi.^2);

g = zeros(size(a));

% All complex numbers -- no information other than a and xi required
g1 = xi + discriminant(~reals);
g2 = xi - discriminant(~reals);
flags = imag(g1)<0;
g1(flags) = g2(flags);
g(~reals) = g1;

% Real numbers -- requires 'branch' information
g(reals) = xi + (2*branch(reals)-1).*discriminant(reals);
