function [eigenvector,eigenvalue] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,N)
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
     
     nSize = 2*N+1; %size
     I=eye(nSize,nSize);
     O=zeros(nSize,nSize);

     % Differentiation Matrices 
     [z1, D1] = fourdif(nSize, 1);
     [z, D2] = fourdif(nSize, 2);
     
     % interval transformation
     z = 0 + (1/(2*pi))*(1-0)*z;

     DX = 1i*kx*I;
     DY = 1i*ky*I;
     DZ = D1;
     Laplacian = -(kx^2+ky^2)*I + D2;

     M1 = diag(-sin(2.*pi*z))*DX + (Pr/Pe)*Laplacian;
     M2 = diag(-sin(2.*pi*z))*DX + (1./Pe)*Laplacian;
     M3 = diag(-sin(2.*pi*z))*DX + (tau/Pe)*Laplacian;
     DU = diag(-2.*pi*cos(2.*pi*z))*I;
     G = (4.*pi*pi*Ri/(Rp-1.))*I;

     A = [I, O, O, O, O, O;
          O, I, O, O, O, O;
          O, O, I, O, O, O;
          O, O, O, O, O, O;
          O, O, O, O, I, O;
          O, O, O, O, O, I];
     B = [M1, O,DU,-DX, O, O;
          O, M1, O,-DY, O, O;
          O, O, M1,-DZ, G,-G;
          DX,DY,DZ, O, O, O;
          O, O, I,  O, M2, O;
          O, O,Rp*I,O, O, M3];

     [eigenvector,eigenvalue] = eig(A,B);
end