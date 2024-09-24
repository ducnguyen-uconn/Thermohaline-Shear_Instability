function MGR_map = maxGrowthRateMap(Ri_list,Pe_list,Rp,Pr,tau)
    sizeOfRi = size(Ri_list,1);
    sizeOfPe = size(Pe_list,1);
    leng = sizeOfRi*sizeOfPe;
    % create an array to store values of growth rate
    MGR = zeros(1,sizeOfRi*sizeOfPe);
    if canUseParallelPool
        % Parallel Computing Toolbox is installed
        parfor index=0:(leng-1)
            Pe_index=floor(index/sizeOfRi);
            Ri_index=mod(index,Pe_index*sizeOfRi);
            % find maximal growth rate
            MGR(index+1) = findMaxGrowthRate(Ri_list(Ri_index+1),Pe_list(Pe_index+1),Rp,Pr,tau);
        end
    else
        % Parallel Computing Toolbox is not installed
        for index=0:(leng-1)
            Pe_index=floor(index/sizeOfRi);
            Ri_index=mod(index,Pe_index*sizeOfRi);
            % find maximal growth rate
            MGR(index+1) = findMaxGrowthRate(Ri_list(Ri_index+1),Pe_list(Pe_index+1),Rp,Pr,tau);
        end
    end
    MGR_map = reshape(MGR,sizeOfRi,sizeOfPe);
end