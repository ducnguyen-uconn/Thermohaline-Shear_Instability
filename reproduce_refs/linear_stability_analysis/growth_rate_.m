clear all;
close all;
clc;


% define paramesters [figure 7b]
Ri_list = [1., 10.]; % Richardson number
Pe_list = [1e4, 1e2]; % Peclet number

Rp = 2.; % density ratio
Pr = 10.; % Prandtl number
tau = 0.01; % diffusivity ratio

k_list=linspace(-0.5,0.5,200);
l_list=linspace(-1,1,200);


N = 128;
I=eye(N,N);
O=zeros(N,N);

% Differentiation Matrices 
[z, D1] = fourdif(N, 1);
[z2, D2] = fourdif(N, 2);

zmin = 0.;
zmax = 1.;
z = zmin+(1./(2.*pi))*(zmax-zmin)*z;

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
                    
                    % compute A
                    A = [I, O, O, O, O, O;
                         O, I, O, O, O, O;
                         O, O, I, O, O, O;
                         O, O, O, O, O, O;
                         O, O, O, O, I, O;
                         O, O, O, O, O, I];

                    % compute B
                    DX = 1i*k*I;
                    DY = 1i*l*I;
                    DZ = D1;
                    Laplacian = -(k^2+l^2)*I + D2;
                    M = diag(-sin(2*pi*z))*DX + (Pr/Pe)*Laplacian;
                    N = diag(-sin(2*pi*z))*DX + (1./Pe)*Laplacian;
                    K = diag(-sin(2*pi*z))*DX + (tau/Pe)*Laplacian;
                    DU = diag(-2*pi*cos(2*pi*z))*I;
                    G = (4.*pi*pi*Ri/(Rp-1.))*I;

                    B = [M, O,DU,-DX, O, O;
                         O, M, O,-DY, O, O;
                         O, O, M,-DZ, G,-G;
                        DX,DY,DZ,  O, O, O;
                         O, O, I,  O, N, O;
                         O, O,Rp*I,O, O, K];
                    
                    % compute eigenvalues
                    [eig_vec, eig_val] = eig(A,B);

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