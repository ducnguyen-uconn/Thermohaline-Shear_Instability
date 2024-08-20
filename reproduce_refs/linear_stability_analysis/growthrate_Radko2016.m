function GR = growthrate_Radko2016(Ri,Pe,Rp,Pr,tau,kx_list,ky_list,N)
     % create an array to store values of growth rate
     GR = zeros(length(kx_list),length(ky_list));
     for kx_index=1:length(kx_list)
          kx=kx_list(kx_index);
          for ky_index=1:length(ky_list)
               ky=ky_list(ky_index);

               % compute eigenvalues
               [eig_vec,eig_val] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,N);
               eig_val(find(real(eig_val)>10^5))=-Inf;
               % compute growth rate
               GR(kx_index,ky_index)=max(real(diag(eig_val)));
          end
     end
end