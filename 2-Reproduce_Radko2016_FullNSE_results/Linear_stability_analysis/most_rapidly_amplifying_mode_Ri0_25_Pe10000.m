clear all;
close all;
clc;
setfigure;

% define paramesters
Ri = 0.25; % Richardson number
Pe = 1e4; % Peclet number
Rp = 2.; % density ratio
Pr = 10.;
tau = 0.01; % diffusivity ratio
kx = 2.68;
ky = 0;
N = 100;%size = 2N+1

[eig_vec,eig_val] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,N);

eig_val(find(real(eig_val)>10^5))=-Inf;

[maxreal,maxindex] = max(real(diag(eig_val)));

varnum = 6; % number of variables [u,v,w,p,T,S]

nSize = 2*N+1;
x=linspace(0.,2.*pi/kx,nSize); % wavelength = 2*pi/kx
z=linspace(0.,1.,nSize);

% create an array to store values of temperature
T = zeros(nSize,nSize);
S = zeros(nSize,nSize);
u = zeros(nSize,nSize);
w = zeros(nSize,nSize);
Tvec = eig_vec(4*nSize+1:4*nSize+nSize,maxindex);
Svec = eig_vec(5*nSize+1:5*nSize+nSize,maxindex);
uvec = eig_vec(0*nSize+1:0*nSize+nSize,maxindex);
wvec = eig_vec(2*nSize+1:2*nSize+nSize,maxindex);

for x_index=1:nSize
    for z_index=1:nSize
        x_val = x(x_index);
        z_val = z(z_index);
        T(x_index,z_index) = real(exp(1i*kx*x_val)*Tvec(z_index));
        S(x_index,z_index) = real(exp(1i*kx*x_val)*Svec(z_index));
        u(x_index,z_index) = real(exp(1i*kx*x_val)*uvec(z_index));
        w(x_index,z_index) = real(exp(1i*kx*x_val)*wvec(z_index));
    end
end
maxT = max(T,[],"all");

% plot 
f1 = figure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcolor(x,z,u'./maxT); shading interp;
colormap(jet);
colorbar;
title('$\hat{u}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk}'})
yticks([0 1])
pbaspect([1 1 1])
savefigure(gca,['most_rapidly_amplifying_mode_Ri=' num2str(Ri) '_Pe=' num2str(Pe) '_N=' num2str(N) '_u.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure;
pcolor(x,z,w'./maxT); shading interp;
colormap(jet);
colorbar;
title('$\hat{w}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk}'})
yticks([0 1])
pbaspect([1 1 1])
savefigure(gca,['most_rapidly_amplifying_mode_Ri=' num2str(Ri) '_Pe=' num2str(Pe) '_N=' num2str(N) '_w.png']);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = figure;
pcolor(x,z,T'./maxT); shading interp;
colormap(jet);
colorbar;
title('$\hat{T}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk}'})
yticks([0 1])
pbaspect([1 1 1])
savefigure(gca,['most_rapidly_amplifying_mode_Ri=' num2str(Ri) '_Pe=' num2str(Pe) '_N=' num2str(N) '_T.png']);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = figure;
pcolor(x,z,S'./maxT); shading interp;
colormap(jet);
colorbar;
title('$\hat{S}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk}'})
yticks([0 1])
pbaspect([1 1 1])
savefigure(gca,['most_rapidly_amplifying_mode_Ri=' num2str(Ri) '_Pe=' num2str(Pe) '_N=' num2str(N) '_S.png']);