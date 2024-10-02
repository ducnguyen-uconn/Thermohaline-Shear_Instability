% define paramesters
Ri = 1.; % Richardson number
Pe = 1e4; % Peclet number
Rp = 2.; % density ratio
Pr = 10.;
tau = 0.01; % diffusivity ratio
kx=2.75;
ky=0;
N = [20 30 50 100 200 300]';

for n_id = 1:length(N)
    n = N(n_id);
    [~,eig_val] = eig_reducedModel(Ri,Pe,Rp,Pr,tau,kx,ky,n);
    eig_val(find(real(eig_val)>1))=-Inf;
    eigv_arr = diag(eig_val)';
    real_matrix = real(eigv_arr);
    [maxvalue,ind] = max(real_matrix,[],"all");
        
    fprintf("N = %d --> max growth rate = %.15f + %.15fi\n",n,real(eigv_arr(ind)),imag(eigv_arr(ind)));
end

% Ri = 1.0, Pe = 1e2, and kx=0.27; 
% N = 50 --> max growth rate = 0.047551258505400 + -0.000000000000028i
% N = 100 --> max growth rate = 0.047551258505366 + 0.000000000000052i
% N = 150 --> max growth rate = 0.047551258505600 + 0.000000000000212i
% N = 200 --> max growth rate = 0.047551258504715 + 0.000000000000470i
% N = 250 --> max growth rate = 0.047551258505086 + 0.000000000000041i
% N = 300 --> max growth rate = 0.047551258502935 + 0.000000000000347i
% N = 350 --> max growth rate = 0.047551258502554 + 0.000000000000905i

% Ri = 10.0, Pe = 1e2, and kx=0.11; 
% N = 50 --> max growth rate = 0.028552449858135 + 0.000000000000062i
% N = 100 --> max growth rate = 0.028552449858150 + -0.000000000000154i
% N = 150 --> max growth rate = 0.028552449858210 + 0.000000000000223i
% N = 200 --> max growth rate = 0.028552449858270 + 0.000000000000297i
% N = 250 --> max growth rate = 0.028552449858906 + 0.000000000000169i
% N = 300 --> max growth rate = 0.028552449858461 + -0.000000000000110i
% N = 350 --> max growth rate = 0.028552449857500 + -0.000000000000056i

% N = 20 --> max growth rate = 0.028552449898713 + 0.000000000000756i
% N = 30 --> max growth rate = 0.028552449858087 + 0.000000000000042i
% N = 50 --> max growth rate = 0.028552449858135 + 0.000000000000062i
% N = 100 --> max growth rate = 0.028552449858150 + -0.000000000000154i
% N = 200 --> max growth rate = 0.028552449858270 + 0.000000000000297i
% N = 300 --> max growth rate = 0.028552449858461 + -0.000000000000110i