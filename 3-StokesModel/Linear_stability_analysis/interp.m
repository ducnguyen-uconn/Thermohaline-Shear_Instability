function [xp,yp,output] = interp(x,y,input,nxp,nyp)
    % x,y,z,xp,yp,zp
    % input - array of input data
    % nxp,nyp,nzp - number of points of new data after interpolating
    % output = zeros(nxp,nyp,nzp);
    xp = linspace(min(x),max(x),nxp); % default: linear
    yp = linspace(min(y),max(y),nyp);

    F = griddedInterpolant({x,y},input,'cubic');
    
    output = F({xp,yp});
end