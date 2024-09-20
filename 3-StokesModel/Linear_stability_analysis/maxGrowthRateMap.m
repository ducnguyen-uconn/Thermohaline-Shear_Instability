% function MGR_map = maxGrowthRateMap(Ri_list,Pe_list,Rp,Pr,tau,N)
%     sizeOfRi = size(Ri_list,1);
%     sizeOfPe = size(Pe_list,1);
%     % create an array to store values of growth rate
%     MGR_map = zeros(sizeOfRi,sizeOfPe);
%     if canUseParallelPool
%         % Parallel Computing Toolbox is installed
%         parfor Ri_index=1:sizeOfRi
%             Ri=Ri_list(Ri_index);
%             for Pe_index=1:sizeOfPe
%                 Pe=Pe_list(Pe_index);
%                 MGR_map(Ri_index,Pe_index) = findMaxGrowthRate(Ri,Pe,Rp,Pr,tau,N);
%             end
%         end
%     else
%         % Parallel Computing Toolbox is not installed
%         for Ri_index=1:sizeOfRi
%             Ri=Ri_list(Ri_index);
%             for Pe_index=1:sizeOfPe
%                 Pe=Pe_list(Pe_index);
%                 MGR_map(Ri_index,Pe_index) = findMaxGrowthRate(Ri,Pe,Rp,Pr,tau,N);
%             end
%         end
%     end
% end


function MGR_map = maxGrowthRateMap(Ri_list,Pe_list,Rp,Pr,tau)
    sizeOfRi = size(Ri_list,1);
    sizeOfPe = size(Pe_list,1);
    leng = sizeOfRi*sizeOfPe;
    % create an array to store values of growth rate
    % MGR_map = zeros(sizeOfRi,sizeOfPe);
    MGR = zeros(1,sizeOfRi*sizeOfPe);
    if canUseParallelPool
        % Parallel Computing Toolbox is installed
        parfor index=0:(leng-1)
            Pe_index=floor(index/sizeOfRi);
            Ri_index=mod(index,Pe_index*sizeOfRi);
            MGR(index+1) = findMaxGrowthRate(Ri_list(Ri_index+1),Pe_list(Pe_index+1),Rp,Pr,tau);
        end
    else
        % Parallel Computing Toolbox is not installed
        for index=0:(leng-1)
            Pe_index=floor(index/sizeOfRi);
            Ri_index=mod(index,Pe_index*sizeOfRi);
            MGR(index+1) = findMaxGrowthRate(Ri_list(Ri_index+1),Pe_list(Pe_index+1),Rp,Pr,tau);
        end
    end
    MGR_map = reshape(MGR,sizeOfRi,sizeOfPe);
end