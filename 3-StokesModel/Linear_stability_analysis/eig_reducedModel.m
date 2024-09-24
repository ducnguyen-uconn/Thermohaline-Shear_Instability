function [eigenvector,eigenvalue] = eig_reducedModel(Ri,Pe,Rp,Pr,tau,kx,ky,N)
     % This function computes eigenvalues of equations in Radko2016's paper 
     % "Thermohaline layering in dynamically and diffusively stable shear flows"
     %
     % Define input paramesters:
     % Ri - Richardson number
     % Pe - Peclet number
     % Rp - density ratio
     % Pr - Prandtl number
     % tau - diffusivity ratio
     % kx, ky - horizontal wavenumbers
     % Define output paramesters:
     % eigenvector and eigenvalues

     % N = 30; % please read the file 'choose_N.m' to see convergence test of N
     
     I=eye(N,N);
     O=zeros(N,N);

     % Differentiation Matrices 
     [~, D1] = fourdif(N, 1);
     [z, D2] = fourdif(N, 2);

     % interval transformation
     z = z/(2*pi);
     D1 = 2*pi*D1;
     D2 = 4*pi*pi*D2;

     DX = 1i*kx*I;
     DY = 1i*ky*I;
     DZ = D1;
     Laplacian = -(kx^2+ky^2)*I + D2;

     M1 = (Pr/Pe)*Laplacian;
     M2 = diag(-sin(2*pi*z))*DX + (1./Pe)*Laplacian;
     M3 = diag(-sin(2*pi*z))*DX + (tau/Pe)*Laplacian;
     G = (4*pi*pi*Ri/(Rp-1.))*I;

     A = [M1, O, O,-DX, O, O;%u
          O, M1, O,-DY, O, O;%v
          O, O, M1,-DZ, G,-G;%w
          DX,DY,DZ, O, O, O;%p
          O, O, I,  O, M2, O;%T
          O, O,Rp*I,O, O, M3];%S
     B = [O, O, O, O, O, O;%u
          O, O, O, O, O, O;%v
          O, O, O, O, O, O;%w
          O, O, O, O, O, O;%p
          O, O, O, O, I, O;%T
          O, O, O, O, O, I];%S

     [eigenvector,eigenvalue] = eig(A,B);
end