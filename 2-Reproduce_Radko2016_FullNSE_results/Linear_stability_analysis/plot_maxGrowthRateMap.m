clear all;
close all;
clc;
setfigure;

% define paramesters [figure 7b]
Ri = logspace(log10(0.25),3,100)'; % Richardson number
Pe = logspace(0,4,100)'; % Peclet number
Rp = 1.5; % density ratio
tau = 0.01; % diffusivity ratio
Pr = 10.;
% filename = ['./test=' num2str(Rp),'_tau=' num2str(tau)];
filename = 'maxGrowthRateMap_RiPe_Rp=2_tau=0.01';

MGR_map = readmatrix([filename '.csv']);
[Rip,Pep,interpMGR_map] = interp_log10(Ri,Pe,MGR_map,100,100);
% vtkwrite([filename '.vtk'],'structured_points','MGR_map', transpose(interpMGR_map));% save to vtk-format file
interpMGR_map(find(interpMGR_map<1e-5))=NaN;
% plot growth rate for each case
% f = figure;
pcolor(Rip,Pep,log10(transpose(interpMGR_map))); colormap(jet);
clim([-6 0])
colorbar;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
xticks([0.25 2 10 200])
yticks([1 10 100 1000 10000])
pbaspect([1 1 1])
title("$\log_{10}$Re$(\lambda)$","Interpreter","latex")
xlabel('{\it{Ri}}')
ylabel('{\it{Pe}}',"Rotation",0)
% savefigure(gca,[filename '.png']);