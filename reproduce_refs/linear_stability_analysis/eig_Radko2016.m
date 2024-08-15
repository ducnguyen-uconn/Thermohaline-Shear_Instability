function [eigenvector,eigenvalue] = eig_Radko2016(Ri,Pe,Rp,tau,kx,ky)
     % This function computes eigenvalues of equations in Radko2016's paper 
     % "Thermohaline layering in dynamically and diffusively stable shear flows"
     %
     % Define input paramesters:
     % Ri - Richardson number
     % Pe - Peclet number
     % Rp - density ratio
     % tau - diffusivity ratio
     % kx, ky - horizontal wavenumbers
     % Define output paramesters:
     % eigenvalues - output eigenvalues
     % 
     % In this function, other parameters are used as
     Pr = 10.; % Prandtl number
     
     N = 128;
     I=eye(N,N);
     O=zeros(N,N);

     % Differentiation Matrices 
     [z1, D1] = fourdif(N, 1);
     [z, D2] = fourdif(N, 2);

     DX = 1i*kx*I;
     DY = 1i*ky*I;
     DZ = D1;
     Laplacian = -(kx^2+ky^2)*I + D2;

     M = diag(-sin(2*pi*z))*DX + (Pr/Pe)*Laplacian;
     N = diag(-sin(2*pi*z))*DX + (1./Pe)*Laplacian;
     K = diag(-sin(2*pi*z))*DX + (tau/Pe)*Laplacian;
     DU = diag(-2*pi*cos(2*pi*z))*I;
     G = (4.*pi*pi*Ri / (Rp-1.))*I;

     A = [I, O, O, O, O, O;
          O, I, O, O, O, O;
          O, O, I, O, O, O;
          O, O, O, O, O, O;
          O, O, O, O, I, O;
          O, O, O, O, O, I];
     B = [M, O,DU,-DX, O, O;
          O, M, O,-DY, O, O;
          O, O, M,-DZ, G,-G;
     DX,DY,DZ,  O, O, O;
          O, O, I,  O, N, O;
          O, O,Rp*I,O, O, K];

     [eigenvector,eigenvalue] = eig(A,B);
end