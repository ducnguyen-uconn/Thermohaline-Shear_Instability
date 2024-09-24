function [xp,yp,output] = interp_log10(x,y,input,nxp,nyp)
    % x,y,z,xp,yp,zp
    % input - array of input data
    % nxp,nyp,nzp - number of points of new data after interpolating
    % output = zeros(nxp,nyp,nzp);
    xp = logspace(log10(min(x)),log10(max(x)),nxp);
    yp = logspace(log10(min(y)),log10(max(y)),nyp);

    F = griddedInterpolant({x,y},input);
    
    output = F({xp,yp});
end