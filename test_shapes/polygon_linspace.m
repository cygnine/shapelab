function z = polygon_range(v,varargin)
% z = polygon_range(v, {step_size=false, N=100})
%
%     Given vertices of a 2D polygon in the complex-valued vector z, returns
%     either: 
%       - a collection of points spaced approximately step_size distance across on
%         each side
%       - a collection of points with exactly N-1 points between vertices
%
%     The "N behavior' overrides the 'step_size' behavior if both are given.
%     Explicity set N to false if you want to use step_size. The polygon is
%     formed in same order as the points are given.

global handles;
opt = handles.common.input_schema({'step_size', 'N'}, {false, 100}, [], varargin{:});

v = v(:);
Nv = length(v);

if opt.N
  z = zeros([opt.N*length(v),1]);
  for q = 1:length(v);
    v1 = v(q);
    v2 = v(mod(q,Nv)+1);
    i1 = (q-1)*opt.N + 1;
    i2 = q*opt.N;

    temp = linspace(v1,v2,opt.N+1).';

    z(i1:i2) = temp(1:end-1);
  end
elseif opt.step_size
  z = [];
  for q = 1:length(v);
    v1 = v(q);
    v2 = v(mod(q,Nv)+1);

    d = norm(v2-v1);
    temp = 0:(opt.step_size)/d:1;
    if temp(end)==1
      temp(end) = [];
    end
    N = length(temp);
    z(end+1:end+N) = temp.'.*(v2-v1) + v1;
  end
else
  error('You seem to have given bad values for both step_size and N');
end
