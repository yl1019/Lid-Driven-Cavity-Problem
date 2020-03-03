clc;
clear;

n = 10;

file1 = fopen('Output.txt','r');
v = fscanf(file1, '%f', [n-2 Inf])';
v_all = zeros(size(v) + 2);

v_all(2:end-1, 2:end-1) = v;
figure(1), contour(rot90(fliplr(v))), axis('square'); 
colorbar;
