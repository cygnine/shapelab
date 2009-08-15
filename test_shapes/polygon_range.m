function Z = polygon_range(V, step_size,ppad)
% z = polygon_range(v, step_size,ppad)
%
%     Given vertices of a 2D polygon in the complex-valued vector V, returns a
%     collection of points spaced approximately step_size distance across on
%     each side. The complex-valued vertices must be specified sequentially.
%
%     Don't know what ppad is about...seems to be amount of forced spacing
%     around vertices. Why would you want this?
%
%     Shamelessly stolen from M. Feiszli.

%ppad = h;
V = V(:);
N = length(V);
Z = [];

if(step_size == 0)
    Z = V;
else
    W = [V;V(1)];

    for n = 1:N
        N = floor(abs(W(n)-W(n+1)) / step_size) + 1;
        tgt  = (W(n+1) - W(n)) / abs(W(n+1) - W(n));
        if(ppad == 0)
            temp = linspace(W(n), W(n+1) - step_size*tgt, N);
        else
            temp = linspace(W(n) + ppad*tgt, W(n+1) - ppad*tgt, N);
        end
        Z = [Z temp];
    end
end
