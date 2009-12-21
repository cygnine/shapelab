function[w] = delta_angle_denoise(z,varargin)
% delta_angle_denoise -- interpolation via delta angle minimzation
%
% w = delta_angle_denoise(z,{k=3,interpolation_weight=10,scale=dyadic})
%
%     Uses linear approximations of beta angles (i.e. the delta angles) to
%     compute interpolations that 'minimize' the norm of the angles. The
%     discrete norm is minimized over k scales, where the definition of 'scale'
%     is determined by the optional input scale. The denoising is done by a
%     least-squares solution since the delta angles are linear functions of the
%     sample points.
%
%     The optional input interpolation_weight specifies how important to
%     consider exactly interpolating the points when solving the least-squares
%     system.
%
%     Adapted from M. Feiszli's code.

persistent indexing delta_angles delta_angle_matrix strict_inputs polygon_linspace
if isempty(strict_inputs)
  from labtools import strict_inputs
  from shapelab.curves import delta_angles delta_angle_matrix indexing
  from shapelab.test_shapes import polygon_linspace
end

inputs = {'k', 'interpolation_weight', 'scale'};
defaults = {3, 10, 'dyadic'};
opt = strict_inputs(inputs, defaults, [], varargin{:});
ind = indexing(opt.scale);

% Parameter determination
z = z(:);
N = length(z);
K = opt.k;
ik = ind(K);

NK = N*ik;

% The matrix we'll try to invert:
D = spalloc(N+NK*K,NK,3*(N+NK*K));

% Coarse scale 
temp = 2*ind(K)*delta_angle_matrix(NK,'k',K,'scale',opt.scale);
D(1:N,:) = opt.interpolation_weight*temp(1:ik:end,:);

% Finer scales
for k = (K-1):-1:0
  temp = 1/2*delta_angle_matrix(NK,'k',k,'scale',opt.scale);
  D((N+(K-k-1)*NK+1):(N+(K-k)*NK),:) = temp;
end

% Restrict preservation of calculated delta angles on coarse scale
b = zeros([N+NK*K 1]);
b(1:N) = ind(1)*opt.interpolation_weight*delta_angles(z, 'k', 0, 'scale', opt.scale);

% Initial guess: linear interpolation
w0 = polygon_linspace(z,'N',ik);

% [x,flag,relres,iter,resvec,lsvec] = lsqr(A,b,tol,maxit,M1,M2,x0)
[w,flag,relres,iter,resvec,lsvec] = lsqr(D,b,1e-8,1000,[],[],w0);
