clear all;
close all;
clc;
setfigure;

% define paramesters [figure 7b]
Ri = 10.; % Richardson number
Pe = 1e2; % Peclet number
Rp = 2.; % density ratio
tau = 0.01; % diffusivity ratio
Pr = 10.;
filename = ['./growth_rate_Ri=' num2str(Ri),'_Pe=' num2str(Pe)];
k_list=linspace(-0.5,0.5,100)';
l_list=linspace(-0.5,0.5,100)';
N = 40; %size=2N+1

GR = growthrate_Radko2016(Ri,Pe,Rp,Pr,tau,k_list,l_list,N);
WriteToVTK(GR, [filename '.vtk']);
GR(find(GR<0))=NaN;
% plot growth rate for each case
f = figure;
pcolor(k_list,l_list,transpose(GR));
% surf(Xq,Zq,transpose(interpGR));view(2);
shading interp; 
colormap(turbo);
colorbar;
pbaspect([1 1 1])
% title(['growth rate Ri=' num2str(Ri),' Pe=' num2str(Pe)])
xlabel('{\it{k}}')
ylabel('{\it{l}}',"Rotation",0)
savefigure(gca,[filename '.png']);