function[] = plot(self)
% plot -- Plots the CircleCurve object
%
% [] = plot(self)

theta = linspace(0, 2*pi, 100);
plot(self.x0(1), self.x0(2), 'x'); hold on;

ang = exp(i*theta);
cx0 = self.x0(1) + i*self.x0(2);

plot(cx0 + self.r*ang)
