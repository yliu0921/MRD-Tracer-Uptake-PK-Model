clc;
clear;close all
% 是不是这样的？

r = (0:.01:1)';

theta = pi*(-1:.01:1);

X = r*cos(theta);

Y = r*sin(theta);

C = X.^2+Y.^2;

p = pcolor(X,Y,C);

set(p,'LineStyle','none');

axis([-1.2 1.2 -1.2 1.2]);

set(gca,'XTick',[],'YTick',[])

axis square