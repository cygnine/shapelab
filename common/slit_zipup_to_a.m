function[z] = slit_zipup_to_a(w,a,varargin)
% [z] = slit_zipup_to_a(w,a,{point_id=zeros(size(z))})
%
%     Implements the inverse of slit_unzip_from_a. This function is defined by
%     Marshall in [1]. Unlike its inverse, this function has an explicit form. 
%
%     0: The point is somewhere in \mathbb{H}\backslash{\mathbb{R}}. Action of
%        the map is normal.
%     1: The point is on \mathbb{R}: care is taken to ensure the output lies on
%        the right line segment
%     2: The point is on the \mathbb{R} between [p-1,p], where p = arg(a)/pi.
%        Although this is an easily-testable case of 1, you may want to force
%        this behavior in the case of machine-epsilon crap.
%
% [1]: Marshall and Rohde, "Convergence of the Zipper algorithm for conformal
%      mapping", 2006.

global handles;
opt = handles.common.InputSchema({'point_id'}, {zeros(size(w))}, [], varargin{:});

interior = opt.point_id==0;
rline = opt.point_id==1;
gamma = opt.point_id==2;

% pre-processing
p = angle(a)/pi;
q = 1-p;
if p<=0 | p>=1
  error('The input a must be in the upper half plane');
end
C = abs(a)/(p^p*q^q);
z = zeros(size(w));

% The 'interior' case
z(interior) = C*(w(interior)-p).^p.*(w(interior)+q).^q;

% The problem cases
z1 = real(w(rline));
z2 = real(w(gamma));
% Let's do some massaging. 'gamma' means on the pre-image of (0,0) -- a
% Find out where the images get zipped up
z1gamma = (-q<=z1) & (z1<=p); % cases where gamma behavior is not explicit (i.e. point_id ~= 2)
z1(z1gamma) = C*(z1(z1gamma)+q).^q.*(p - z1(z1gamma)).^p.*exp(i*pi*p);

% On other locations of the real line:
temp = z1(~z1gamma);
z1(~z1gamma) = C*sign(temp).*abs(temp-p).^p.*abs(temp+q).^q;
z(rline) = z1;

% If things are on gamma, straightforward
z2 = C*(z2+q).^q.*(p-z2).^p.*exp(i*pi*p);
z(gamma) = z2;
