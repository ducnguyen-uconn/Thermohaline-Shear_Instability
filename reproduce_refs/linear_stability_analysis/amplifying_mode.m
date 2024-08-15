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
kx = 0.125;
ky = 0;
N = 1;%size = 2N+1

[eig_vec,eig_val] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,N);

fprintf('size of eig_vec is %s\n', mat2str(size(eig_vec)));

[maxreal,maxindex] = max(real(diag(eig_val)));

x=linspace(0,2*pi/kx,101); % wavelength = 2*pi/kx
z=linspace(0,1,101);

% create an array to store values of temperature
T = zeros(length(x),length(z));
for x_index=1:length(x)
    for z_index=1:length(z)
        exp1 = exp(1i*kx*x(x_index) + 0);
        sumexp2 = sum(eig_vec,"all");
        T(x_index,z_index) = real(exp1*sumexp2);
    end
end

% plot 
f = figure;
pcolor(x,z,T); shading interp;
colormap(turbo);
colorbar;
title('$\hat{T}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk}'})
yticks([0 1])
pbaspect([2 1 1])
savefigure(gca,'amplifying_mode_test.png');