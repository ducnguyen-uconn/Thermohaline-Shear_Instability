% define paramesters
Ri = 0.25; % Richardson number
Pe = 1e4; % Peclet number
Rp = 2.; % density ratio
Pr = 10.;
tau = 0.01; % diffusivity ratio
kx=5.66;
ky=0;
N = [50 100 150 200 250 300 400 500]';

for n_id = 1:length(N)
    n = N(n_id);
    [~,eig_val] = eig_Radko2016(Ri,Pe,Rp,Pr,tau,kx,ky,n);
    eig_val(find(real(eig_val)>1))=-Inf;
    eigv_arr = diag(eig_val)';
    real_matrix = real(eigv_arr);
    [maxvalue,ind] = max(real_matrix,[],"all");
        
    fprintf("N = %d --> max growth rate = %.15f + %.15fi\n",n,real(eigv_arr(ind)),imag(eigv_arr(ind)));
end

% Ri = 1.0, Pe = 1e4, and kx=2.75; represents a case of KH instability
% N = 30 --> max growth rate = 0.056048087224452 + -0.285830106315358i
% N = 50 --> max growth rate = 0.028757442473861 + 0.172127191093144i
% N = 60 --> max growth rate = 0.000672583477145 + 6.243487777354595i
% N = 70 --> max growth rate = 0.010237044521566 + 0.123090599065128i
% N = 80 --> max growth rate = 0.000672583477143 + 6.243487777354598i
% N = 90 --> max growth rate = 0.000672583477187 + -6.243487777354583i
% N = 100 --> max growth rate = 0.000672583477167 + -6.243487777354606i  <----
% N = 150 --> max growth rate = 0.000672583477126 + 6.243487777354614i
% N = 200 --> max growth rate = 0.000672583477119 + 6.243487777354574i
% N = 300 --> max growth rate = 0.000672583477081 + -6.243487777354466i

% Ri = 10.0, Pe = 1e2, and kx=0.11; represents a case without KH instability ---> rapid grid convergence
% N = 30 --> max growth rate = 0.028033138652798 + 0.000000000000016i
% N = 50 --> max growth rate = 0.028033138652797 + -0.000000000000076i
% N = 60 --> max growth rate = 0.028033138652870 + -0.000000000000020i
% N = 70 --> max growth rate = 0.028033138652781 + 0.000000000000104i
% N = 80 --> max growth rate = 0.028033138652451 + -0.000000000000136i
% N = 90 --> max growth rate = 0.028033138652354 + -0.000000000000416i
% N = 100 --> max growth rate = 0.028033138653061 + 0.000000000000069i
% N = 150 --> max growth rate = 0.028033138653730 + 0.000000000000183i
% N = 200 --> max growth rate = 0.028033138653368 + 0.000000000001967i

% Ri = 0.25, Pe = 1e4, and kx=5.66; represents a case with strong KH instability ---> convergence beginning 250
% N = 50 --> max growth rate = 0.246237614132212 + -0.271687125447877i
% N = 100 --> max growth rate = 0.123568713323971 + -0.140065302333401i
% N = 150 --> max growth rate = 0.051182637276454 + 0.091877801796247i
% N = 200 --> max growth rate = -0.018731641322890 + 0.060650025089352i
% N = 250 --> max growth rate = -0.033871040380110 + 5.636969829659659i
% N = 300 --> max growth rate = -0.033871040380114 + 5.636969829659674i
% N = 400 --> max growth rate = -0.033871040380065 + 5.636969829659686i
% N = 500 --> max growth rate = -0.033871040380046 + -5.636969829659686i