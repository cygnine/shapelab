function[U] = kernel_basis(theta, varargin)
% kernel_basis -- Evaluates elements in ker(L)
%
% U = kernel_basis(theta, [[d=0]])
%
%     This function standardizes the basis description of ker(L). It takes in a
%     vector theta, and the matrix U has three columns: each contains
%     evaluations of a basis function at the locations theta. Currently, the
%     following basis is employed:
%
%     U = [1, sin, cos]
%
%     The input d determines how many derivatives of the basis functions to
%     take.

if nargin<2
  d = 0;
else
  d = round(varargin{1});
end

theta = theta(:);
N = length(theta);
U = zeros([N 3]);
if d>0
  U(:,1) = 0;
else
  U(:,1) = ones(size(theta));
end

if mod(d,2)==0
  U(:,2) = sin(theta);
  U(:,3) = cos(theta);
  U(:,2:3) = U(:,2:3)*(-1)^(d/2);
else
  U(:,2) = cos(theta);
  U(:,3) = -sin(theta);
  U(:,2:3) = U(:,2:3)*(-1)^((d-1)/2);
end
