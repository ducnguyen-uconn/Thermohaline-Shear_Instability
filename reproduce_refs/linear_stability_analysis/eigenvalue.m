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

N = 64;
I=eye(N,N);
O=zeros(N,N);

% Differentiation Matrices 
[z, D1] = fourdif(N, 1);
[z2, D2] = fourdif(N, 2);

DX = 1i*kx*I;
DY = 1i*ky*I;
DZ = D1;
Laplacian = -(kx^2+ky^2)*I + D2;

% change to diag()*D
M = -sin(2*pi*z).*DX + (Pr/Pe)*Laplacian;
N = -sin(2*pi*z).*DX + (1./Pe)*Laplacian;
K = -sin(2*pi*z).*DX + (tau/Pe)*Laplacian;
DU = -2*pi*cos(2*pi*z).*I;
G = 4.*pi*pi*Ri / (Rp-1.)*I;

A = [I, O, O, O, O, O;
     O, I, O, O, O, O;
     O, O, I, O, O, O;
     O, O, O, O, O, O;
     O, O, O, O, I, O;
     O, O, O, O, O, I];
B = [M, O,DU,-DX, O, O;
     O, M, O,-DY, O, O;
     O, O, M,-DZ, G,-G;
    DX,DY,DZ,  O, O, O;
     O, O, I,  O, N, O;
     O, O,Rp*I,O, O, K];

[eig_vec, eig_val] = eig(A,B);


% plot eigenvalues
f = figure;
plot(real(eig_val),imag(eig_val),'o');
line([0 0], ylim,'Color','black');  %x-axis
line(xlim, [0 0],'Color','black');  %y-axis
xlabel('real') 
ylabel('imag') 
saveas(f,'./eig_val_Real','png')