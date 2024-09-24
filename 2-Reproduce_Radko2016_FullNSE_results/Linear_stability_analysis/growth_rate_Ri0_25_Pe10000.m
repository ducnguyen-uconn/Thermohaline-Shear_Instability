clear all;
close all;
clc;
setfigure;

% define paramesters [figure 7b]
Ri = 0.25; % Richardson number
Pe = 1e4; % Peclet number
Rp = 2.; % density ratio
tau = 0.01; % diffusivity ratio
Pr = 10.;
filename = ['./growth_rate_Ri=' num2str(Ri),'_Pe=' num2str(Pe)];
k_list=linspace(-10,10,200)';
l_list=linspace(-5,5,100)';
N=200;
GR = growthrate_Radko2016(Ri,Pe,Rp,Pr,tau,k_list,l_list,N);

[kp,lp,interpGR] = interp(k_list,l_list,GR,1000,1000);
[maxmax,idx] = max(interpGR,[],"all");
[kmax,lmax] = ind2sub(size(interpGR),idx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <---- export data in vtk-format
% vtkwrite([filename '.vtk'],'structured_points','growth_rate', transpose(interpGR));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interpGR(find(interpGR<0))=NaN; % NaN values will be displayed on figure -> white color
% plot growth rate
f = figure;
pcolor(kp,lp,transpose(interpGR));shading interp;colormap(jet);colorbar;
pbaspect([1 0.5 1])
title(['max=' num2str(maxmax) 'at kx=' num2str(kmax) ',kx=' num2str(kmax)])
xlabel('{\it{k}}')
ylabel('{\it{l}}',"Rotation",0)
savefigure(gca,[filename '.png']);