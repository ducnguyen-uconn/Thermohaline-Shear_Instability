clear all;
close all;
clc;

% define paramesters [figure 7]
Ri_list = [1.]; % Richardson number
Pe_list = [1e4]; % Peclet number
Rp = 2.; % density ratio
Pr = 10.; % Prandtl number
tau = 0.01; % diffusivity ratio

k_list=linspace(-4,4,201);
l_list=logspace(-4,4,401);
phi_list=1;

z=1;%????????????????????????????????

for Ri_index=1:length(Ri_list)
     Ri=Ri_list(Ri_index);
     for Pe_index=1:length(Pe_list)
          Pe=Pe_list(Pe_index);
          for k_index=1:length(k_list)
               k=k_list(k_index);
               for l_index=1:length(l_list)
                    l=l_list(l_index);
                    
                    % compute A
                    A = [1, 0, 0, 0, 0, 0;
                         0, 1, 0, 0, 0, 0;
                         0, 0, 1, 0, 0, 0;
                         0, 0, 0, 0, 0, 0;
                         0, 0, 0, 0, 1, 0;
                         0, 0, 0, 0, 0, 1];

                    % compute B
                    DX = 1i*k;
                    DY = 1i*l;
                    DZ = 0;%?
                    Laplacian = -(k^2+l^2) + 0;%?
                    M = -sin(2*pi*z)*DX + (Pr/Pe)*Laplacian;%<---------
                    N = -sin(2*pi*z)*DX + (1./Pe)*Laplacian;
                    K = -sin(2*pi*z)*DX + (tau/Pe)*Laplacian;
                    DU = -2*pi*cos(2*pi*z);
                    G = 4.*pi*pi*Ri / (Rp-1.);

                    B = [M, 0,DU,-DX, 0, 0;
                         0, M, 0,-DY, 0, 0;
                         0, 0, M,-DZ, G,-G;
                        DX,DY,DZ,  0, 0, 0;
                         0, 0, 1,  0, N, 0;
                         0, 0,Rp,  0, 0, K];
                    
                    % compute eigenvalues
                    [eig_vec, eig_val] = eig(A,B);

                    %compute growth rate
                    growth_rate{Ri_index,Pe_index}(k_index,l_index)=max(real(diag(eig_val)));
               end
          end
     end
end

growth_rate{1}(find(growth_rate{1}<0))=NaN;
data{1}.x=k_list;
data{1}.y=l_list;
data{1}.z=growth_rate{1}';
plot_config.print_size=[1,1000,900];
plot_config.name=['growth_rate_Ri=',num2str(Ri_list),'_Pe=',num2str(Pe_list),'.png'];
plot_config.label_list={1,'$k$','$l$'};
plot_contour(data,plot_config);
