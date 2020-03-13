clf;clear;clc;

data = dlmread('Output.txt', '', 1, 0);
data = sortrows(data);
x = data(:, 1);
y = data(:, 2);
v = data(:, 3);
s = data(:, 4);

dist    = abs(y - 0.5);
minDist = min(dist);
idx     = find(dist == minDist);
vx = data(idx, 5);

dist    = abs(x - 0.5);
minDist = min(dist);
idy     = find(dist == minDist);
vy = data(idy, 6);


x = unique(x); nx = length(x);
y = unique(y); ny = length(y);
[X, Y] = meshgrid(x, y);
vorticity = reshape(v, ny, nx);
streamfunction = reshape(s, ny, nx);

% plot contour
figure(1)
subplot(1, 2, 1);
contourf(X, Y, vorticity);
title('vorticity'); axis('square');
subplot(1, 2, 2);
contourf(X, Y, streamfunction);
title('stream function'); axis('square');

% plot velocity
% figure(2)
% plot(x, vx);
% figure(3)
% plot(y, vy);