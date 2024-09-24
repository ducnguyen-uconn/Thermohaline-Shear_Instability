clear all;
close all;
clc;
setfigure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change this
Ri = 0.25; % Richardson number
Pe = 2; % Peclet number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change resolution here
N_list=[40,60,80,100]';
legend_list = "N="+num2str(N_list);
filename = ['./convergence_test_of_max_growth_rate_Ri=' num2str(Ri),'_Pe=' num2str(Pe)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed parameters
Rp = 2.; % density ratio
tau = 0.01; % diffusivity ratio
Pr = 10.;
k_list=linspace(0.00001,10,50)';
l_list=0;
f = figure;
for N_index=1:length(N_list)
    N = N_list(N_index);
    fprintf('N=%d\n',N);
    GR = growthrate_reducedModel(Ri,Pe,Rp,Pr,tau,k_list,l_list,N);
    % [maxmax,kmax] = max(GR,[],"all");
    plot(k_list,GR','Color', [0 0 0 N_index/length(N_list)],'LineWidth',2);hold on;
    % plot(k_list(kmax),maxmax,'o','Color', [0 0 0 N_index/length(N_list)]);hold on;
end
legend(legend_list);
% pbaspect([1 0.5 1])
xlabel('{\it{k}}')
ylabel('Growth rate')
savefigure(gca,[filename '.png']);