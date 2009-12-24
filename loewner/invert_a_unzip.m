function [gplus,gminus] = invert_a_unzip(aplus,aminus, xi);
% invert_a_unzip -- inverts the g(a) relation with intention for unzipping
%
% [gplus,gminus] = invert_a_unzip(aplus,aminus, xi)
%
%     The relation a = g^2 - 2*g*xi is inverted to solve for g, with appropriate
%     choices of branches for g. The branch in the upper-half plane is always
%     chosen.
%
%     When a is complex, the branch is always chosen automatically; when a is
%     real, both branches are computed.
%
%     When aplus is complex, aminus should be as well. They should only differ
%     when they are both real-valued.
%
%     The two outputs gplus and gminus correspond to both branches of the
%     inversion relation. When a is complex, the values are the same; when a is
%     real-valued, different branches are chosen.

reals = (imag(aplus)==0);

discriminant = sqrt(2*aplus + xi.^2);

gplus = zeros(size(aplus));

% All complex numbers -- no information other than a and xi required
g1 = xi + discriminant(~reals);
g2 = xi - discriminant(~reals);
flags = imag(g1)<0;
g1(flags) = g2(flags);
gplus(~reals) = g1;
gminus = gplus;

% Real numbers -- requires 'branch' information
gplus(reals) = xi + discriminant(reals);
gminus(reals) = xi - sqrt(2*aminus(reals) + xi.^2);
