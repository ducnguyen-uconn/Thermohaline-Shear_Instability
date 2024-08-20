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

eig_val(find(real(eig_val)>10^12))=-Inf;

% writematrix(eig_vec,'eig_vec_tab.txt','Delimiter','tab')
% writematrix(eig_val,'eig_val_tab.txt','Delimiter','tab')
% fprintf('size of eig_vec is %s\n', mat2str(size(eig_vec)));

[maxreal,maxindex] = max(real(diag(eig_val)))

varnum = 6; % number of variables [u,v,w,p,T,S]


x=linspace(0,2*pi/kx,101); % wavelength = 2*pi/kx
z=linspace(0,1,101);

% create an array to store values of temperature
T = zeros(length(x),length(z));
S = zeros(length(x),length(z));
u = zeros(length(x),length(z));
w = zeros(length(x),length(z));
for x_index=1:length(x)
    for z_index=1:length(z)
        exp1 = exp(1i*kx*x(x_index) + 0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sumexp2 = 0;
        for index=-N:N
            sumexp2 = sumexp2+eig_vec(maxindex,(index+N)*varnum+1)*exp(2*pi*1i*index*z(z_index));
        end
        u(x_index,z_index) = real(exp1*sumexp2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sumexp2 = 0;
        % for index=-N:N
        %     sumexp2 = sumexp2+eig_vec(maxindex,(index+N)*varnum+3)*exp(2*pi*1i*index*z(z_index));
        % end
        % w(x_index,z_index) = real(exp1*sumexp2);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sumexp2 = 0;
        % for index=-N:N
        %     sumexp2 = sumexp2+eig_vec(maxindex,(index+N)*varnum+5)*exp(2*pi*1i*index*z(z_index));
        % end
        % T(x_index,z_index) = real(exp1*sumexp2);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sumexp2 = 0;
        % for index=-N:N
        %     sumexp2 = sumexp2+eig_vec(maxindex,(index+N)*varnum+6)*exp(2*pi*1i*index*z(z_index));
        % end
        % S(x_index,z_index) = real(exp1*sumexp2);
    end
end

% plot 
f1 = figure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcolor(x,z,u); shading interp;
% colormap(turbo);
colorbar;
title('$\hat{u}$','Interpreter','latex')
xlabel('{\it{x}}')
ylabel('{\it{z}}',"Rotation",0)
xticks([0 2*pi/kx])
xticklabels({'0','2\pi/{\itk}'})
yticks([0 1])
pbaspect([1 1 1])
savefigure(gca,'amplifying_mode_test_u.png');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f2 = figure;
% pcolor(x,z,w); shading interp;
% % colormap(turbo);
% colorbar;
% title('$\hat{w}$','Interpreter','latex')
% xlabel('{\it{x}}')
% ylabel('{\it{z}}',"Rotation",0)
% xticks([0 2*pi/kx])
% xticklabels({'0','2\pi/{\itk}'})
% yticks([0 1])
% pbaspect([2 1 1])
% savefigure(gca,'amplifying_mode_test_w.png');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f3 = figure;
% pcolor(x,z,T); shading interp;
% % colormap(turbo);
% colorbar;
% title('$\hat{T}$','Interpreter','latex')
% xlabel('{\it{x}}')
% ylabel('{\it{z}}',"Rotation",0)
% xticks([0 2*pi/kx])
% xticklabels({'0','2\pi/{\itk}'})
% yticks([0 1])
% pbaspect([2 1 1])
% savefigure(gca,'amplifying_mode_test_T.png');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f4 = figure;
% pcolor(x,z,S); shading interp;
% % colormap(turbo);
% colorbar;
% title('$\hat{S}$','Interpreter','latex')
% xlabel('{\it{x}}')
% ylabel('{\it{z}}',"Rotation",0)
% xticks([0 2*pi/kx])
% xticklabels({'0','2\pi/{\itk}'})
% yticks([0 1])
% pbaspect([2 1 1])
% savefigure(gca,'amplifying_mode_test_S.png');