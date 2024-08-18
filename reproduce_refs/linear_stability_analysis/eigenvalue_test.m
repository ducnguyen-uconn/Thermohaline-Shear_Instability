clear all;
close all;
clc;
setfigure;

% define paramesters
Ri = 10.; % Richardson number
Pe = 1e2; % Peclet number
Rp = 2.; % density ratio
Pr = 10.;
tau = 0.01; % diffusivity ratio
kx = 0.1;
ky = 0;
N = 100;%size = 2N+1

[eig_vec,eig_val] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,N);

eig_val=diag(eig_val);
fprintf('max of Re(eig_val) is %s\n', mat2str(max(real(eig_val))));

% plot eigenvalues
f = figure;
plot(real(eig_val),imag(eig_val),'o');
line([0 0], ylim,'Color','black');  %x-axis
line(xlim, [0 0],'Color','black');  %y-axis
xlabel('real') 
ylabel('imag',"Rotation",0) 
title(['eigenvalues Ri=' num2str(Ri),' Pe=' num2str(Pe)])
savefigure(gca,['eigenvalues_Ri=' num2str(Ri) '_Pe=' num2str(Pe) '_N=' num2str(N) '.png']);