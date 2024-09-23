clear all;
close all;
clc;
setfigure;

% define paramesters [figure 7b]
Ri = 1.; % Richardson number
Pe = 1e4; % Peclet number
Rp = 2.; % density ratio
tau = 0.01; % diffusivity ratio
Pr = 10.;
filename = ['./growth_rate_Ri=' num2str(Ri),'_Pe=' num2str(Pe)];
k_list=linspace(-4,4,100)';
l_list=linspace(-4,4,100)';

GR = growthrate_Radko2016(Ri,Pe,Rp,Pr,tau,k_list,l_list,80);
[kp,lp,interpGR] = interp(k_list,l_list,GR,1000,1000);
[maxmax,idx] = max(interpGR,[],"all");
[kmax,lmax] = ind2sub(size(interpGR),idx);
disp(['max value = ' num2str(maxmax) ' of growth rate at (k,l)=(' num2str(abs(kp(kmax))) ',' num2str(abs(lp(lmax))) ')']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% <---- export data in vtk-format
% vtkwrite([filename '.vtk'],'structured_points','growth_rate', transpose(interpGR));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interpGR(find(interpGR<0))=NaN; % NaN values will be displayed on figure -> white color
% plot growth rate
f = figure;
pcolor(kp,lp,transpose(interpGR));shading interp;colormap(jet);colorbar;
pbaspect([1 1 1])
xlabel('{\it{k}}')
ylabel('{\it{l}}',"Rotation",0)
savefigure(gca,[filename '.png']);