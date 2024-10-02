function maxGR = findMaxGrowthRate(Ri,Pe,Rp,Pr,tau)
    % adaptive resolutions based on Pe
    optN = 40;
    if 1e2<=Pe && Pe<1e3
        optN = 100;
    elseif 1e3<=Pe
        optN = 200;
    end

    % adaptive wavenumber ranges based on Ri
    opt_max_kx = 4.0;
    if 1<=Ri && Ri<10
        opt_max_kx = 0.5;
    elseif 10<=Ri
        opt_max_kx = 0.15;
    end
    kx_list=linspace(0.0001,opt_max_kx,50)';
    sizeOfkx = size(kx_list,1);
    ky=0;

    maxGR = -inf;
    for kx_index=1:sizeOfkx
        kx=kx_list(kx_index);
        % compute eigenvalues
        [~,eig_val] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,optN);
        eig_val(find(real(eig_val)>10^5))=-Inf;
        % compute growth rate
        if max(real(diag(eig_val))) > maxGR
            maxGR = max(real(diag(eig_val)));
        end
    end
end