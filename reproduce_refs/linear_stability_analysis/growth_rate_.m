clear all;
close all;
clc;
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'defaultTextFontName', 'Times New Roman');

% define paramesters [figure 7b]
Ri_list = [1., 10.]; % Richardson number
Pe_list = [1e4, 1e2]; % Peclet number

Rp = 2.; % density ratio
tau = 0.01; % diffusivity ratio

k_list=linspace(-0.5,0.5,10);
l_list=linspace(-1,1,10);


for Ri_index=1:length(Ri_list)
     Ri=Ri_list(Ri_index);
     for Pe_index=1:length(Pe_list)
          Pe=Pe_list(Pe_index);

          % create an array to store values of growth rate
          GR = zeros(length(k_list),length(l_list));

          for k_index=1:length(k_list)
               k=k_list(k_index);
               for l_index=1:length(l_list)
                    l=l_list(l_index);
                    
                    % compute eigenvalues
                    [eig_vec,eig_val] = eig_Radko2016(Ri,Pe,Rp,tau,k,l);

                    % compute growth rate
                    GR(k_index,l_index)=max(real(diag(eig_val)));
               end
          end

          % plot growth rate for each case
          f = figure;
          % pcolor(k_list,l_list,GR); 
          % pcolor(k_list,l_list,GR); shading flat;
          pcolor(k_list,l_list,GR); shading interp;
          colorbar
          savingname = ['./growth_rate_Ri=' num2str(Ri),'_Pe=' num2str(Pe)];
          title(['growth rate Ri=' num2str(Ri),' Pe=' num2str(Pe)])
          saveas(f,savingname,'png')

          % WriteToVTK(GR, [savingname '.vtk']);
     end
end