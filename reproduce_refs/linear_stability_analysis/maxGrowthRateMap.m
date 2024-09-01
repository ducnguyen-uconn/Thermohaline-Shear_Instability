function MGR_map = maxGrowthRateMap(Ri_list,Pe_list,Rp,Pr,tau,N)
    kx_list=linspace(0.01,4,50)';
    ky_list=0;
    sizeOfkx = size(kx_list,1);
    sizeOfky = size(ky_list,1);
    sizeOfRi = size(Ri_list,1);
    sizeOfPe = size(Pe_list,1);
    % create an array to store values of growth rate
    MGR_map = zeros(sizeOfRi,sizeOfPe);
    if canUseParallelPool
        % Parallel Computing Toolbox is installed
        parfor Ri_index=1:sizeOfRi
            Ri=Ri_list(Ri_index);
            for Pe_index=1:sizeOfPe
                Pe=Pe_list(Pe_index);
                maxGR = -inf;
                for kx_index=1:sizeOfkx
                    kx=kx_list(kx_index);
                    ky=ky_list;
                    % compute eigenvalues
                    [eig_vec,eig_val] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,N);
                    eig_val(find(real(eig_val)>10^5))=-Inf;
                    % compute growth rate
                    if max(real(diag(eig_val))) > maxGR
                        maxGR = max(real(diag(eig_val)));
                    % else
                    %     break;
                    end
                end
                MGR_map(Ri_index,Pe_index) = maxGR;
            end
        end
    else
        % Parallel Computing Toolbox is not installed
        for Ri_index=1:sizeOfRi
            Ri=Ri_list(Ri_index);
            for Pe_index=1:sizeOfPe
                Pe=Pe_list(Pe_index);
                maxGR = -inf;
                for kx_index=1:sizeOfkx
                    kx=kx_list(kx_index);
                    ky=ky_list;
                    % compute eigenvalues
                    [eig_vec,eig_val] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,N);
                    eig_val(find(real(eig_val)>10^5))=-Inf;
                    % compute growth rate
                    if max(real(diag(eig_val))) > maxGR
                        maxGR = max(real(diag(eig_val)));
                    %else
                    %    break;
                    end
                end
                MGR_map(Ri_index,Pe_index) = maxGR;
            end
        end
    end
end