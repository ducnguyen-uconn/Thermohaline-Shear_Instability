clear all;
close all;
clc;

% define paramesters
Ri = 1.; % Richardson number
Rp = 2.; % density ratio
Pe = 1e4; % Peclet number
Pr = 10.; % Prandtl number
tau = 0.01; % diffusivity ratio

kx = 2;
ky = 0;

N = 8;
I=eye(N,N);
O=zeros(N,N);

% Differentiation Matrices 
[z, D1] = fourdif(N, 1);
[z2, D2] = fourdif(N, 2);

DX = kx*I;
DY = ky*I;
DZ = D1;
Laplacian = -(kx^2+ky^2)*I + D2;

M = -sin(2*pi*z).*DX + (Pr/Pe)*Laplacian;
N = -sin(2*pi*z).*DX + (1./Pe)*Laplacian;
K = -sin(2*pi*z).*DX + (tau/Pe)*Laplacian;
DU = -2*pi*cos(2*pi*z);
G = 4.*pi*pi*Ri / (Rp-1.)*I;

A = diag([I I I O I I])
B = [M, O,DU,-DX, O, O;
     O, M, O,-DY, O, O;
     O, O, M,-DZ, G,-G;
    DX,DY,DZ,  O, O, O;
     O, O, I,  O, N, O;
     O, O,Rp,  O, O, K]

[eig_vec, eig_val] = eig(A,B);
