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

[maxreal,maxindex] = max(real(diag(eig_val)));

x=linspace(0,2*pi/kx,101);
z=linspace(0,1,101);

% create an array to store values of growth rate
T = zeros(length(x),length(z));
for x_index=1:length(x)
    for z_index=1:length(z)
        T(x_index,z_index) = real(exp(1i*kx*x(x_index) + 1i*z(z_index))*sum(eig_vec*exp(2*1i*pi*z),"all"));
    end
end


% plot 
f = figure;
pcolor(x,z,T); shading interp;
colorbar
saveas(f,'./amplifying_mode_test','png')