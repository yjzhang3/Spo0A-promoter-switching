%% this file is trying to find how RNAp (more specificlaly, sigma-H bound
% holonenzyme concentration increaeses with time 
load('promoter_activity_single.mat')
load('Spo0A_dynamics.mat')
real_data_Ps = Ps_promoter_activity_mean(2:7,1:8); % exclude the real WT, instead WT is Ps only, no mutations

TF_conc_t = wtaCopy(1:6:end);
RNAp_conc_t_mock = zeros(length(TF_conc_t),1)+1;

nbd = 5;
mut_mat = [[1,1,1];[0,1,1];[1,0,1];[1,1,0];[0,0,1];[0,1,0];[0,0,1];[0,0,0]];
energy_len = nbd+nbd*(nbd-1)/2;

%% bounds

lb = zeros(energy_len,1)-20;
ub = zeros(energy_len,1)+15;

lb(1:4) = 5; % standard energy must be positive

% binding box 2 and 4 have lower affinity (high standard energy)
lb(2) = 9;
lb(4) = 6.5; 

lb(5) = 4; % promoter energy should be no greater than standard

 % introduce TF-TF repulsions
lb([6:8 10:11 13]) = 1; 
lb(11) = 5; % box 2 and 4 repulsion is larger
ub([6:8 10:11 13]) = 10; % TF-repulsion shouldn't be dominate

lb(9) = 5; % make one of them necessarily a repressor
ub(9) = 10; % but not a super strong repressor
%% fit 123*  data
% treat both eenergies and RNAp as unknowns

for ll = 1:length(real_data_Ps(1,:))
    real_data_Ps(:,ll) = real_data_Ps(:,ll)./max(real_data_Ps(:,ll));
end

lb_overall = zeros(energy_len+6,1);
ub_overall = zeros(energy_len+6,1);
lb_overall(1:energy_len) = lb;
ub_overall(1:energy_len) = ub;

lb_overall(energy_len+1:end) = zeros(length(TF_conc_t),1);
ub_overall(energy_len+1:end) = zeros(length(TF_conc_t),1)+1000000;

all_mut_data = real_data_Ps(:,end); % the 123* strain
mut_RNAp = [0,0,0];

%% find [RNAp] first
[new_pars,diff] = fit_data_sigma(nbd,TF_conc_t,mut_RNAp,all_mut_data,lb_overall,ub_overall);

% %% then find energy
% RNAp_conc_t_real = new_pars(energy_len+1:end);
% energyi_sigma = new_pars(1:energy_len);
% [final_energyi,diff] = fit_data_new(nbd,TF_conc_t,RNAp_conc_t_real,mut_mat,real_data_Ps,lb,ub);

%% now see  how the newe parameters do
RNAp_conc_t_real = new_pars(energy_len+1:end);
sigma_energyi = new_pars(1:energy_len);

figure();
for kk = 1:length(real_data_Ps)
    subplot(4,2,kk)
     
    plot(TF_conc_t,real_data_Ps(:,kk),'o','LineWidth',4)
    ylim([0 1])
    hold on
    
    TR = time_dep_TR_new_wSigma(nbd,sigma_energyi,TF_conc_t,RNAp_conc_t_real,mut_mat(kk,:));
    plot(TF_conc_t,TR,'LineWidth',4)
    xlabel('TF concentration')
    ylabel('transcription rate')
    title(sprintf('mutation %d',kk))
end