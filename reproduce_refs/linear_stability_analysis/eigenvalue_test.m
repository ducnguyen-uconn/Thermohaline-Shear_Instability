clear all;
close all;
clc;
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'defaultTextFontName', 'Times New Roman');

% define paramesters
Ri = 1.; % Richardson number
Pe = 1e4; % Peclet number
Rp = 2.; % density ratio
tau = 0.01; % diffusivity ratio
kx = 2*pi/2;
ky = 0;

[eig_vec,eig_val] = eig_Radko2016(Ri,Pe,Rp,tau,kx,ky);

max(real(diag(eig_val)))

% plot eigenvalues
f = figure;
plot(real(eig_val),imag(eig_val),'o');
line([0 0], ylim,'Color','black');  %x-axis
line(xlim, [0 0],'Color','black');  %y-axis
xlabel('real') 
ylabel('imag') 
saveas(f,'./eigenvalues_test','png')