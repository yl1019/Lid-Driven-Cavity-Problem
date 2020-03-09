clf;clear;clc;

data = dlmread('Output.txt', '', 1, 0);
data = sortrows(data);
x = data(:, 1);
y = data(:, 2);
v = data(:, 3);
s = data(:, 4);



x = unique(x); nx = length(x);
y = unique(y); ny = length(y);
[X, Y] = meshgrid(x, y);
vorticity = reshape(v, ny, nx);
streamfunction = reshape(s, ny, nx);

subplot(1, 2, 1);
contourf(X, Y, vorticity);
title('vorticity'); axis('square');
subplot(1, 2, 2);
contourf(X, Y, streamfunction);
title('stream function'); axis('square');