function[psi] = weld_primitive_extend_driver(x, stuff)
% weld_primitive_extend_driver -- Evaluates the real-line extension of the weld primitive
%
% psi = weld_primitive_extend_driver(x, stuff)
%
%     Uses the information from the struct stuff to evaluate the primitive of
%     the weld. stuff is the output from
%     weld_primitive_setup. More details can be found in
%     weld_primitive_extend.

persistent weval
if isempty(weval)
  from shapelab.extensions import weld_primitive_evaluate as weval
end

psi = zeros(size(x));
psi2pi = stuff.global_primitive; % saves typing

xnegative = x<0;
xpositive = not(xnegative);

% For x positive it's not so bad:
shifts = floor(x(xpositive)/(2*pi));

psi(xpositive) = psi2pi*shifts + (4*pi^2)*mysum(shifts) + ...
                 shifts*2*pi.*mod(x(xpositive), 2*pi) + ...
                 weval(mod(x(xpositive), 2*pi), stuff);

% For x negative we need to be a bit more creative:
%shifts = ceil(x(xnegative)/(2*pi));
shifts = floor(x(xnegative)/(2*pi));
% "I0 of negative 2 pi"
I0n2pi = 4*pi^2 - psi2pi;

xtemp = mod(x(xnegative), 2*pi);

psi(xnegative) = 2*pi*(2*pi - xtemp) - (psi2pi - weval(xtemp, stuff)) - ...
                 I0n2pi*(shifts+1) - (4*pi^2)*mysum(shifts) - ...
                 (shifts+1)*2*pi.*(2*pi - xtemp);
end

function[result] = mysum(n)
% mysum -- a silly subfunction
  result = zeros(size(n));
  gt1 = n>1;
  result(gt1) = 1/2*(n(gt1)-1).*(abs(n(gt1)-1)+1);

  lt2 = n<-2;
  result(lt2) = 1/2*(n(lt2)+2).*(abs(n(lt2)+2) +1);
end
