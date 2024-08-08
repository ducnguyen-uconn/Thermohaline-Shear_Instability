clear all;
close all;
clc;

% define paramesters
Ri = 1.; % Richardson number
Rp = 2.; % density ratio
Pe = 1e4; % Peclet number
Pr = 10.; % Prandtl number
tau = 0.01; % diffusivity ratio

N = 64;

[x, D] = fourdif(N, 1);
