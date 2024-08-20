clear all;
close all;
clc;
setfigure;

% define paramesters [figure 7b]
Ri_list = [10.]; % Richardson number
Pe_list = [1e2]; % Peclet number

Rp = 2.; % density ratio
tau = 0.01; % diffusivity ratio

k_list=linspace(-0.5,0.5,50);
l_list=linspace(-0.5,0.5,50);

Pr = 10.;

N = 30; %size=2N+1


for Ri_index=1:length(Ri_list)
     Ri=Ri_list(Ri_index);
     for Pe_index=1:length(Pe_list)
          Pe=Pe_list(Pe_index);

          GR = growthrate_Radko2016(Ri,Pe,Rp,Pr,tau,k_list,l_list,N);
          GR(find(GR<0))=0;
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % [X,Z] = meshgrid(k_list,l_list);
          % [Xq,Zq] = meshgrid(-0.5:0.01:0.5);
          % interpGR = interp2(X,Z,GR,Xq,Zq,'cubic');
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % plot growth rate for each case
          f = figure;
          pcolor(k_list,l_list,transpose(GR));
          % surf(Xq,Zq,transpose(interpGR));view(2);
          shading interp; 
          colormap(turbo);
          colorbar;
          pbaspect([1 1 1])
          filename = ['./growth_rate_Ri=' num2str(Ri),'_Pe=' num2str(Pe)];
          title(['growth rate Ri=' num2str(Ri),' Pe=' num2str(Pe)])
          xlabel('{\it{k}}')
          ylabel('{\it{l}}',"Rotation",0)
          savefigure(gca,[filename '.png']);

          % WriteToVTK(GR, [filename '.vtk']);
     end
end